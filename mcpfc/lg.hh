#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <unistd.h>
#include <pthread.h>
#include <sys/time.h>
#include <csignal>
#include <acml.h>

namespace GLX{
#include <GL/glx.h>
}

#include <boost/program_options.hpp>

class Display{
public:
  Display(int*, char**, int, int);
  ~Display();
  void keyPressed(unsigned char, int, int);
  void mouseAction(int, int, int, int);
  void DrawGLScene();
  void startEventLoop();

private:
  static const int m_escape = 27;
  int m_window;
  int m_wWidth;
  int m_wHeight;
  int m_scale;
  unsigned int m_texture[1];
  int m_tWidth, m_tHeight, m_tElements;
  unsigned char* m_tData;
  GLuint m_base;

  void setUpOpenGL();
  void setUpTexture();
  void BuildFont();
  void glPrint(const char*);
};

class Simulation{
public:
  Simulation(int, double, double, bool, double, double, double, double, int);
  ~Simulation();
  double nBorSum(double*, int, int);
  double nBorSum2(double*, int, int);
  double nBorSum3(double*, int, int);
  void run();
  int m_width, m_height;
  double m_initProb;

  double* m_domain;
  double* m_varPsiPrim;
  double* m_varPsiSec;
  double* m_varLapPsi;
  double* m_varLapLapPsi;
  double* m_varEnergy;
  double m_mu;
  double m_r;
  double m_tempK;
  double m_size;
  int m_ms;

  double m_tol;
  int* m_state;
  int m_info;
  bool m_run;
  bool m_novis;

  pthread_mutex_t m_mutexData;
};
