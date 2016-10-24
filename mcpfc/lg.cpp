#include "lg.hh"

using namespace std;
namespace po = boost::program_options;

Display* display;
Simulation* s;
unsigned int ctr = 0;

void DrawGLScene(){ display->DrawGLScene(); }
void keyPressed(unsigned char key, int x, int y){ display->keyPressed(key, x, y); }
void mouseAction(int button, int state, int x, int y){ display->mouseAction(button, state, x, y); }
void Idle(){ display->DrawGLScene(); usleep(10000); }
void* simThread(void* ptr){ s->run(); }

inline int wrap(const int n, const int w){
  if(n < 0) return n + w;
  else if(n >= w) return n - w;
  else return n;
}

void doubleToRGB(double& r, double& g, double& b, double h, double min, double max, double s = 1.0, double v = 1.0){
  if(fabs(max - min) < 1e6){ r = 0.0; g = 0.0; b = 0.0;}
  h-=min; h/=(max - min); h*=360.0;
  double x = s * v * (1.0 - fabs(fmod(h / 60.0, 2.0) - 1.0));
  switch((int)(h / 60.0)){
  case 0: r = s * v; g = x; b = 0.0; break;
  case 1: r = x; g = s * v; b = 0.0; break;
  case 2: r = 0.0; g = s * v; b = x; break;
  case 3: r = 0.0; g = x; b = s * v; break;
  case 4: r = x; g = 0.0; b = s * v; break;
  default: r = s * v; g = 0.0; b = x;
  }
}

void Display::BuildFont(){
  m_base = glGenLists(96);
  GLX::Display* dpy = GLX::XOpenDisplay(NULL);
  GLX::XFontStruct* fontInfo = GLX::XLoadQueryFont(dpy, "-misc-fixed-*-*-*-*-20-*-*-*-*-*-*-*");
  //  GLX::XFontStruct* fontInfo = GLX::XLoadQueryFont(dpy, "-*-helvetica-*-r-*-*-20-*-*-*-*-*-*-*");
  if(fontInfo == NULL){
    fontInfo = GLX::XLoadQueryFont(dpy, "fixed");
    if (fontInfo == NULL) {
      printf("no X font available?\n");
    }
  }
  GLX::glXUseXFont(fontInfo->fid, 32, 96, m_base);
  GLX::XFreeFont(dpy, fontInfo);
  GLX::XCloseDisplay(dpy);
}

void Display::glPrint(const char *text){
  if(text == NULL){
    return;
  }
  glPushAttrib(GL_LIST_BIT);
  glListBase(m_base - 32);
  glCallLists(strlen(text), GL_UNSIGNED_BYTE, text);
  glPopAttrib();
}

void Display::setUpTexture(){
  glEnable( GL_TEXTURE_RECTANGLE_NV );
  glGenTextures(1, &m_texture[0]);
  glBindTexture( GL_TEXTURE_RECTANGLE_NV, m_texture[0]);
  glTexParameteri(GL_TEXTURE_RECTANGLE_NV,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
  glTexParameteri(GL_TEXTURE_RECTANGLE_NV,GL_TEXTURE_MIN_FILTER,GL_NEAREST_MIPMAP_NEAREST);
  // 2d texture, level of detail 0 (normal), 3 components (red, green, blue), x size from image, y size from image,
  // border 0 (normal), rgb color data, unsigned byte data, and finally the data itself.
  glTexImage2D(GL_TEXTURE_RECTANGLE_NV, 0, 4, m_tWidth, m_tHeight, 0, GL_RGBA, GL_UNSIGNED_BYTE, m_tData);
}

void Display::DrawGLScene(){
  pthread_mutex_lock(&s->m_mutexData);
  int len = s->m_height * s->m_width;
  double* dom = s->m_domain;
  double min = 1e6;
  double max = -1e6;
  for(int siteI = 0; siteI < len; siteI++){
    double value = dom[siteI];
    if(value < min)
      min = value;
    else if(value > max)
      max = value;
  }
  //  printf("min: %f, max: %f\n",min,max);
  for(int i=0,j=0; i < len; i++,j+=4){
    double value = dom[i];
    double r, g, b;
    doubleToRGB(r, g, b, value, min, max);
    m_tData[j] = 255 * r;
    m_tData[j + 1] = 255 * g;
    m_tData[j + 2] = 255 * b;
  }
  // for(int ro = 0; ro < m_scale; ro++){
  //   for(int co = 0; co < m_scale; co++){
  //      int tmpIndex = 3 * (startIndex + ro * m_tWidth + co);
  pthread_mutex_unlock(&s->m_mutexData);
  // glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  // glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
  // glLoadIdentity();
  glBindTexture(GL_TEXTURE_RECTANGLE_NV, m_texture[0]);
  glTexSubImage2D(GL_TEXTURE_RECTANGLE_NV, 0, 0, 0, m_tWidth, m_tHeight, GL_RGBA, GL_UNSIGNED_BYTE, m_tData);
  glBegin(GL_QUADS);
  // origo: top-left
  // glTexCoord2f( 0, m_tHeight); glVertex2f( -1, -1);
  // glTexCoord2f( m_tWidth, m_tHeight); glVertex2f( 1, -1);
  // glTexCoord2f( m_tWidth, 0); glVertex2f( 1, 1);
  // glTexCoord2f( 0, 0); glVertex2f( -1, 1);

  // origo: bottom-left
  glTexCoord2i( 0, 0); glVertex2i( -1, -1);
  glTexCoord2i( m_tWidth, 0); glVertex2i( 1, -1);
  glTexCoord2i( m_tWidth, m_tHeight); glVertex2i( 1, 1);
  glTexCoord2i( 0, m_tHeight); glVertex2i( -1, 1);
  glEnd();

  //void* font = GLUT_BITMAP_9_BY_15;
  char text[100];
  sprintf(text, "T=%.4f, min: %.4f, max: %.4f", s->m_tempK, min, max);

  glDisable(GL_TEXTURE_RECTANGLE_NV);
  glPushAttrib(GL_CURRENT_BIT);
  //  glColor3f(0.0,0.0,1.0);
  glColor3f(1.0,1.0,1.0);
  glRasterPos2f(-0.95,-0.95);
  glPrint(text);
  // int i=0;
  // while(text[i]!=0){
  //   glutBitmapCharacter(font, text[i]); i++;
  // }
  glPopAttrib();
  glEnable(GL_TEXTURE_RECTANGLE_NV);

  // string s = "Respect mah authoritah!";
  // void * font = GLUT_BITMAP_9_BY_15;
  // for (string::iterator i = s.begin(); i != s.end(); ++i){
  //   char c = *i;
  //   glutBitmapCharacter(font, c);
  // }


  glutSwapBuffers();
}

void saveToFile(){
  pthread_mutex_lock(&s->m_mutexData);
  char fname[50];
  sprintf(fname, "dom%09d.dat", ctr);
  FILE* out = fopen(fname, "wb");
  if(out == NULL){
    printf("couldnt open file %s for writing\n",fname);
    pthread_mutex_unlock(&s->m_mutexData);
    return;
  }
  int len = s->m_height * s->m_width;
  double* dom = s->m_domain;
  fwrite(dom, sizeof(double), len, out);
  fclose(out);
  printf("file %s saved\n",fname);
  pthread_mutex_unlock(&s->m_mutexData);
}

void Display::keyPressed(unsigned char key, int x, int y){
  //  printf("key: %d\n", key);
  switch(key){
  case m_escape:
    glutDestroyWindow(m_window);
    exit(0);
    break;
  case 's':
    saveToFile();
    break;
  }
  usleep(1000);
}

void Display::mouseAction(int button, int state, int x, int y){
  double& t = s->m_tempK;
  double DT = 0.05 * t;
  if(button == 3 && state == 0){
    t += DT;
  }else if(button == 4 && state == 0){
    t -= DT;
  }
  //  printf("t: %f\n", t);
  //printf("m: %d : %d\n", button, state);
  usleep(1000);
}

void Display::setUpOpenGL(){
  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_ALPHA);
  glutInitWindowSize(m_tWidth, m_tHeight);
  glutInitWindowPosition(0, 0);
  m_window = glutCreateWindow("MC");
  glutDisplayFunc(&::DrawGLScene);
  glutIdleFunc(&Idle);
  glutKeyboardFunc(&::keyPressed);
  glutMouseFunc(&::mouseAction);
  setUpTexture();
  glEnable(GL_TEXTURE_RECTANGLE_NV);
  glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
  gluOrtho2D(-1, 1, -1, 1);
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  //  glPushAttrib(GL_DEPTH_BUFFER_BIT | GL_LIGHTING_BIT);
  glDisable(GL_LIGHTING);
  glDisable(GL_DITHER);
  glDisable(GL_BLEND);
  glDisable(GL_DEPTH_TEST);
  glEnable(GL_COLOR_MATERIAL);
  glViewport(0, 0, m_tWidth, m_tHeight);
  BuildFont();
}

void Display::startEventLoop(){
  glutMainLoop();
}

Display::Display(int* argcp, char** argv, int width, int height){
  m_scale = 1;
  m_tWidth = m_scale * width;
  m_tHeight = m_scale * height;
  m_tElements = m_tWidth * m_tHeight;
  m_tData = new unsigned char [4 * m_tElements];

  glutInit(argcp, argv);
  setUpOpenGL();
}

Display::~Display(){
  delete[] m_tData;
  glDeleteLists(m_base, 96);
}

Simulation::Simulation(int divs, double size, double prob, bool novis, double temp, double tol, double r, double mu, int ms){
  m_ms = ms;
  m_size = size;
  m_r = r;
  m_mu = mu;
  m_tol = tol;
  m_novis = novis;
  m_run = true;
  m_mutexData = PTHREAD_MUTEX_INITIALIZER;
  m_width = divs;
  m_height = divs;
  m_initProb = prob;
  m_tempK = temp;
  srand(time(NULL));

  int elements = m_width * m_height;

  m_varPsiPrim = new double[elements];
  m_varPsiSec = new double[elements];
  m_varLapPsi = new double[elements];
  m_varLapLapPsi = new double[elements];
  m_varEnergy = new double[elements];

  m_domain = m_varPsiPrim;

  for(int i = 0; i < m_width * m_height; i++){
    //    double p = rand() / (double) RAND_MAX;
    m_varPsiPrim[i] = m_initProb; //2.0 * (p - 0.5) * m_initProb;
  }

  int genid = 3;
  int subid = 0;
  int seed[1]; //MT providing 1 instead of 624 required
  seed[0] = rand();
  int lstate = 633; // genid == 3; MT lstate >= 633
  int lseed = 1; //only 1 i supplied, others generated
  m_state = new int[lstate];
  // drandinitialize(genid, subid, seed, &lseed, m_state, &lstate, &m_info);
  if(m_info != 0) printf("prng was not initialized successfully\n");
}

Simulation::~Simulation(){
  delete[] m_varPsiPrim;
  delete[] m_varPsiSec;
  delete[] m_varLapPsi;
  delete[] m_varLapLapPsi;
  delete[] m_varEnergy;

  delete[] m_state;
}

void printTime(timeval startTime, timeval endTime){
  printf("%d", (int)(((endTime.tv_sec  - startTime.tv_sec) * 1000
                      + (endTime.tv_usec - startTime.tv_usec)/1000.0)
                     + 0.5));
}

// double Simulation::nBorSum(double* ar, int r, int c, double* o1, double* o2, double* o3){
//   double res=0;
//   int rp = wrap(r + 1, m_height);
//   int rm = wrap(r - 1, m_height);
//   int cp = wrap(c + 1, m_width);
//   int cm = wrap(c - 1, m_width);
//   res += ar[rp * m_width + c];
//   res += ar[rm * m_width + c];
//   res += ar[r * m_width + cp];
//   res += ar[r * m_width + cm];

//   // res += m_domain[rp * m_width + cp];
//   // res += m_domain[rp * m_width + cm];
//   // res += m_domain[rm * m_width + cp];
//   // res += m_domain[rm * m_width + cm];

//   return res;
// }

// double Simulation::nBorSum2(double* ar, int r, int c){
//   double res=0;
//   int rp = wrap(r + 1, m_height);
//   int rm = wrap(r - 1, m_height);
//   int cp = wrap(c + 1, m_width);
//   int cm = wrap(c - 1, m_width);
//   res += ar[rp * m_width + cp];
//   res += ar[rm * m_width + cm];
//   res += ar[rm * m_width + cp];
//   res += ar[rp * m_width + cm];
//   return res;
// }

// double Simulation::nBorSum3(double* ar, int r, int c){
//   double res=0;
//   int rp = wrap(r + 2, m_height);
//   int rm = wrap(r - 2, m_height);
//   int cp = wrap(c + 2, m_width);
//   int cm = wrap(c - 2, m_width);
//   res += ar[rp * m_width + c];
//   res += ar[rm * m_width + c];
//   res += ar[r * m_width + cp];
//   res += ar[r * m_width + cm];
//   return res;
// }

void Simulation::run(){
  double* pick = new double[m_width * m_height];
  //  double* choose = new double[m_width * m_height];
  double* sign = new double[m_width * m_height];

  unsigned int counter = 0;
  timeval startTime, endTime;
  gettimeofday(&startTime, NULL);

  while(m_run && (m_ms == -1 || ctr <= m_ms)){
    double dh = m_size / m_width;
    double dh2 = dh * dh;
    double dh4 = dh2 * dh2;
    // laplace psi
    // for(int r = 0; r < m_height; r++){
    //   for(int c = 0; c < m_width; c++){
    //     m_varLapPsi[r * m_width + c] = nBorSum(m_varPsiPrim, r, c) - 4.0 * m_varPsiPrim[r * m_width + c];
    //   }
    // }
    // // laplace^2 psi
    // for(int r = 0; r < m_height; r++){
    //   for(int c = 0; c < m_width; c++){
    //     m_varLapLapPsi[r * m_width + c] = nBorSum(m_varLapPsi, r, c) - 4.0 * m_varLapPsi[r * m_width + c];
    //   }
    // }
    // energy
    // for(int r = 0; r < m_height; r++){
    //   for(int c = 0; c < m_width; c++){
    //     m_varEnergy[r * m_width + c] =
    //       double psi = m_varPsiPrim[r * m_width + c];
    //       - m_mu * psi + 0.5 * (m_r + 1.0) * psi * psi + 0.25 * psi * psi * psi * psi +
    //         psi * m_varLapPsi[r * m_width + c] / dh2 + 0.5 * psi * m_varLapLapPsi[r * m_width + c] / (dh2 * dh2);
    //   }
    // }
    // Delta E / Delta Psi
    // for(int r = 0; r < m_height; r++){
    //   for(int c = 0; c < m_width; c++){
    //     double psi = m_varPsiPrim[r * m_width + c];
    //     m_varEnergy[r * m_width + c] =
    //       - m_mu + (m_r + 1.0) * psi + psi * psi * psi +
    //       m_varLapPsi[r * m_width + c] / dh2 + 0.5 * m_varLapLapPsi[r * m_width + c] / (dh2 * dh2);
    //   }
    // }

    //    if(not m_novis){
      // if(counter < 10000000){
      //   counter+= 3 * m_batchSize;
      // }else{
      //   counter = 0;
      //   gettimeofday(&endTime, NULL);
      //   cout << "10M random sampling took: ";
      //   printTime(startTime, endTime);
      //   cout << " ms" << endl;
      //   gettimeofday(&startTime, NULL);
      // }
    ctr++;
    if(ctr % 10 == 0){
      double rho = 0;
      for(int siteI = 0; siteI < m_height * m_width; siteI++){
        rho += m_varPsiPrim[siteI];
      }
      rho /= (m_height * m_width);
      printf("sweep: %d, avg.Psi: %.5f\n",ctr,rho);
    }
    /*
    if(ctr % 100 == 0){
      saveToFile();
    }
    */

    //    }
    // ctr++;
    // if(ctr % 10 == 0){
    //   double magnetization = 0;
    //   double energyPerSite = 0;
    //   for(int siteI = 0; siteI < m_height * m_width; siteI++){
    //     int r = siteI / m_width;
    //     int c = siteI - r * m_width;
    //     magnetization += m_domainPrim[siteI];
    //     energyPerSite += (- m_domainPrim[siteI] * nBorSum(r, c));
    //   }
    //   energyPerSite /= (double)(m_width * m_height);
    //   magnetization /= (double)(m_width * m_height);
    //   if(m_novis){
    //     printf("%d\t%.5f\t%.5f\n",ctr,energyPerSite,magnetization);
    //   }else{
    //     printf("ctr: %d, sites: %d, avg.mag.: %.5f, avg.energy: %.5f\n",ctr, m_height * m_width, magnetization, energyPerSite);
    //   }
    // }

    //dranduniform(m_width * m_height, 0.0, 1.0, m_state, pick, &m_info);
    //    dranduniform(m_width * m_height, 0.0, 1.0, m_state, choose, &m_info);
    //dranduniform(m_width * m_height, 0.0, 1.0, m_state, sign, &m_info);
    for(int i = 0; i < m_width * m_height; i++){
        pick[i] = rand() / (double)RAND_MAX;
        sign[i] = rand() / (double)RAND_MAX;
    }

    // CORRECT
    // double p4 = exp(- 4.0 / m_tempK);
    // double p8 = exp(- 8.0 / m_tempK);
    // CORRECT

    // FUN
    //    int nbors[4][2] = {{+1,0},{-1,0},{0,+1},{0,-1}};
    //    int nbors[8][2] = {{+1,0},{-1,0},{0,+1},{0,-1},{+1,+1},{+1,-1},{-1,-1},{-1,+1}};
    // FUN

    pthread_mutex_lock(&m_mutexData);
    for(int siteI = 0; siteI < m_width * m_height; siteI++){
      int r = siteI / m_width;
      int c = siteI - r * m_width;

      int rp = wrap(r + 1, m_height);
      int rm = wrap(r - 1, m_height);
      int cp = wrap(c + 1, m_width);
      int cm = wrap(c - 1, m_width);

      int rpp = wrap(r + 2, m_height);
      int rmm = wrap(r - 2, m_height);
      int cpp = wrap(c + 2, m_width);
      int cmm = wrap(c - 2, m_width);

      double nsum1=
        m_varPsiPrim[rp * m_width + c] +
        m_varPsiPrim[rm * m_width + c] +
        m_varPsiPrim[r * m_width + cp] +
        m_varPsiPrim[r * m_width + cm];

      double nsum2=
        m_varPsiPrim[rp * m_width + cp] +
        m_varPsiPrim[rp * m_width + cm] +
        m_varPsiPrim[rm * m_width + cp] +
        m_varPsiPrim[rm * m_width + cm];

      double nsum3=
        m_varPsiPrim[rpp * m_width + c] +
        m_varPsiPrim[rmm * m_width + c] +
        m_varPsiPrim[r * m_width + cpp] +
        m_varPsiPrim[r * m_width + cmm];

      double oldPsi = m_varPsiPrim[siteI];
      double newPsi = oldPsi + 2.0 * (sign[siteI] - 0.5) * m_tol;

      double oldPsi2 = oldPsi * oldPsi;
      double newPsi2 = newPsi * newPsi;

      double dE = 0.5 * (m_r + 1.0) * (newPsi2 - oldPsi2) + 0.25 * (newPsi2 * newPsi2 - oldPsi2 * oldPsi2) - m_mu * (newPsi - oldPsi) +
        2.0 / dh2 * (newPsi - oldPsi) * (nsum1 - 2.0 * newPsi - 2.0 * oldPsi) +
        1.0 / dh4 * (newPsi - oldPsi) * (- 8.0 * nsum1 + 2.0 * nsum2 + nsum3 + 10.0 * (newPsi + oldPsi));
      dE *= dh2;

      double newPsiPrime = oldPsi;

      if(dE <= 0.0){
        newPsiPrime = newPsi;
      }else{
        double p = expf(- dE / m_tempK);
        if(pick[siteI] < p){
          newPsiPrime = newPsi;
        }
      }
      m_varPsiSec[siteI] = newPsiPrime;
    }
    // exchange domains
    double* tmpdomain = m_varPsiPrim;
    m_varPsiPrim = m_varPsiSec;
    m_varPsiSec = tmpdomain;

    pthread_mutex_unlock(&m_mutexData);
  }

  delete[] pick;
  //  delete[] choose;
  delete[] sign;
}

void ctrlCHandler(int){
  s->m_run = false;
  cerr << " terminating" << endl;
}

int main(int argc, char** argv){
  // parameters:
  bool novis = false;
  double initprob = 0.5;
  double temp = 2.4;
  int divs = 100;
  double size = 2 * M_PI;
  double tol = 0.2;
  double r = -0.25;
  double mu = 0.0;
  int ms = -1;

  try{
    po::options_description desc("Options");
    desc.add_options()
      ("help", "produce help message")
      ("iprob", po::value<double>(), "initial spin up probability")
      ("novis", "turn off visualization")
      ("divs", po::value<int>(), "use a divs x divs discretization")
      ("size", po::value<double>(), "simulate a size x size domain")
      ("tol", po::value<double>(), "site pick tolerance")
      ("temp", po::value<double>(), "temperature")
      ("mu", po::value<double>(), "mu")
      ("ms", po::value<int>(), "max sweeps")
      ("r", po::value<double>(), "r");
      //      ("prefix", po::value<double>(), "prefix for saved files");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    if (vm.count("help")) { cout << desc << "\n"; return 0; }
    if (vm.count("iprob")){
      double option = vm["iprob"].as<double>();
      initprob = option;
      //      if(option <= 1.0 && option >= 0.0) initprob = option;
    }
    if (vm.count("size")){
      double option = vm["size"].as<double>();
      if(option >= 0.5) size = option;
    }
    if (vm.count("divs")){
      int option = vm["divs"].as<int>();
      if(option >= 10) divs = option;
    }
    if (vm.count("ms")){
      int option = vm["ms"].as<int>();
      if(option > 0) ms = option;
    }
    if (vm.count("temp")){
      double option = vm["temp"].as<double>();
      if(option > 0.0) temp = option;
    }
    if (vm.count("tol")){
      double option = vm["tol"].as<double>();
      if(option >= 0.0 && option <= 1.0) tol = option;
    }
    if (vm.count("r")){
      double option = vm["r"].as<double>();
      if(option >= -1.0 && option <= 0.0) r = option;
    }
    if (vm.count("mu")){
      double option = vm["mu"].as<double>();
      mu = option;
      //if(option >= -2.0 && option <= 2.0) mu = option;
    }
    if (vm.count("novis")) novis = true;
  }
  catch(exception& e) { cerr << "error: " << e.what() << "\n"; return 1; }
  catch(...) { cerr << "Exception of unknown type!\n"; return 1; }

  cerr << "parameters:"
       << "\ntemperature: " << temp
       << "\ninitial spin up probability: " << initprob
       << "\ndomain size: " << size
       << "\nr: " << r
       << "\nmu: " << mu
       << "\nsite pick tolerance: " << tol
       << "\nmax sweeps: " << ms << endl;

  signal(SIGINT, ctrlCHandler);

  s = new Simulation(divs, size, initprob, novis, temp, tol, r, mu, ms);
  if(not novis) display = new Display(&argc, argv, divs, divs);
  pthread_t simTh; pthread_create(&simTh, NULL, simThread, NULL);
  if(not novis) display->startEventLoop();
  pthread_join(simTh, NULL);
  delete s; delete display;
  return 0;
}

// g++ main.cpp -O3 -std=c++0x -I ~/acml/gfortran64/include/ -L ~/acml/gfortran64/lib/ -lX11 -lGLU -lGL -lglut -lm -lpthread -lacml -ffast-math -fexpensive-optimizations -o lg

// for n in `seq 7 20`; do cd ..; mkdir "run"$n; cd "run"$n; ../../lg --divs 100 --r -0.2 --mu -0.2 --size 100 --tol 0.1 --temp 1.0 --iprob 0 --novis --ms 50000; done


// for NN in `seq 17 40`; do cd ..; mkdir "run"$NN; cd "run"$NN; ../../lg --divs 1000 --r -0.2 --mu -0.2 --size 1000 --tol 0.1 --temp 0.4 --iprob 0 --novis --ms 16000; done
