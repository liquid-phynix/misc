#include<iostream>
#include<vector>
#include<list>
#include<map>
#include<set>
#include<algorithm>
#include<eigen3/Eigen/Dense>
#include<cmath>
#include<cstdio>
#include<complex>

using namespace std;
using namespace Eigen;

typedef Vector3d Vec;
typedef Vector3i IVec;
typedef Matrix3d Mat;
typedef Matrix3i IMat;
typedef pair<IVec, double> Dir;

struct DirCompare{
  bool operator()(Dir da, Dir db){
    auto& a=da.first; auto& b=db.first;
    if(a[0] < b[0]) return true;
    else if(a[0] > b[0]) return false;
    else if(a[1] < b[1]) return true;
    else if(a[1] > b[1]) return false;
    else return a[2] < b[2];
  }
};

typedef list<Dir> DirList;
typedef set<Dir, DirCompare> DirSet;

inline bool isAngleLess(IVec iv1, IVec iv2){
  Vec v1=iv1.cast<double>(); Vec v2=iv2.cast<double>();
  return acos(v1.dot(v2)/sqrt(v1.dot(v1)*v2.dot(v2))) < (4.0*M_PI_2/5.0);
}

inline bool isAngleLess(Vec v1, Vec v2){
  return acos(v1.dot(v2)/sqrt(v1.dot(v1)*v2.dot(v2))) < (4.0*M_PI_2/5.0);
}

DirSet genDirections(DirList& dirs){
  int indices[]={0,1,2};
  DirList tmpDirs;
  do{
    for(auto& d: dirs)
      tmpDirs.emplace_back(Dir(IVec((d.first)(indices[0]),(d.first)(indices[1]),(d.first)(indices[2])), d.second));
  }while(next_permutation(indices, indices+3));
  list<IVec> signs={{1,1,1},{-1,1,1},{1,-1,1},{1,1,-1},{-1,-1,1},{-1,1,-1},{1,-1,-1},{-1,-1,-1}};
  DirSet allDirs;
  for(auto& sign: signs)
    for(auto& dir: tmpDirs)
      allDirs.insert(Dir(sign.array() * dir.first.array(), dir.second));
  return allDirs;
}

list<Vec> genTriplesAndSolve(DirSet& allDirs){
  MatrixXi boolmat(allDirs.size(), allDirs.size());
  int row=0, col;
  for(auto it1=allDirs.begin(); it1!=allDirs.end(); it1++, row++){
    auto it2=it1; it2++; col=row+1;
    for(; it2!=allDirs.end(); it2++, col++)
      boolmat(col, row) = boolmat(row, col) = isAngleLess(it1->first, it2->first);
  }

  list<Vec> points; int tcnt=0;
  auto it1=allDirs.begin(); int i1=0;
  for(; it1!=allDirs.end(); it1++, i1++){
    auto it2=it1; it2++; int i2=i1+1;
    for(; it2!=allDirs.end(); it2++, i2++){
      auto it3=it2; it3++; int i3=i2+1;
      for(; it3!=allDirs.end(); it3++, i3++){
        if(boolmat(i1,i2) and boolmat(i1,i3) and boolmat(i2,i3)){
          tcnt++;
          Mat m;
          m.row(0)=it1->first.cast<double>();
          m.row(1)=it2->first.cast<double>();
          m.row(2)=it3->first.cast<double>();
          Vec b(it1->second * m.row(0).norm(),
                it2->second * m.row(1).norm(),
                it3->second * m.row(2).norm());
          Mat im=m.inverse();
          Vec point=im*b;
          //cerr << "norm: " << (m*point-b).norm()/b.norm() << "\n";
          double relErr=(m*point-b).norm()/b.norm();
          double np=point.norm();
          if(relErr<1e-3 and np>3.0 and np<5.0)
            points.push_back(point);
        }
      }
    }
  }
  cerr << "triples: " << tcnt << "\n";
  return points;
}

bool isOutside(Vec tvec, Vec rvec){
  double a2=tvec.dot(tvec);
  Vec tmr=tvec-rvec;
  double b2=tmr.dot(tmr);
  double arg=(a2+b2-rvec.dot(rvec))/(2.0*sqrt(a2)*sqrt(b2));
  return real(acos(complex<double>(arg, 0))) > (M_PI_2+0.001);
}

list<Vec> erasePoints(DirSet& allDirs, list<Vec>& points){
  list<Vec> gammas;
  for(auto& d: allDirs){
    Vec v=d.first.cast<double>();
    v.array() *= d.second / v.norm();
    gammas.push_back(v);
  }
  list<Vec> ret;
  for(auto& p: points){
    bool include=true;
    for(auto& g: gammas){
      if(isAngleLess(g, p) and isOutside(g,p)){
        include=false;
        break;
      }
    }
    if(include) ret.push_back(p);
  }
  return ret;
}

int main(){
  DirList data;
  while(not cin.eof()){
    int in[3]; double g;
    cin >> in[0] >> in[1] >> in[2] >> g >> ws;
    sort(in, in+3);
    data.push_back(Dir(IVec(in[0], in[1], in[2]), g));
  }
  
  auto allDirs=genDirections(data);
  cerr << "allDirs done" << endl; cerr << "allDirs length: " << allDirs.size() << endl;
  
  auto points=genTriplesAndSolve(allDirs);
  cerr << "triples and points done" << endl; cerr << "points length: " << points.size() << endl;
  
  auto newPoints=erasePoints(allDirs, points);
  cerr << "erasing done" << endl; cerr << "remained points: " << newPoints.size() << endl;
   
  cout << 3 << "\n" << newPoints.size() << "\n";
  for(auto it=newPoints.begin(); it!=newPoints.end(); it++){
    printf("%.10f %.10f %.10f\n", (*it)(0), (*it)(1), (*it)(2));
    //cout << it->t();
    //cout << (*it)(0) << " " <<  (*it)(1) << " " << (*it)(2) << "\n";
  }
  cout << endl;

  return 0;
}

// g++ -std=c++0x wulff.cpp -O3 -o wulffcpp
// cat gamma.original | ./wulffcpp > points.from.cpp
// cat points.from.cpp | qhull QJ m > polygons.from.cpp
// cat points.from.cpp | qhull C0.01 Qx C-0 m > polygons.from.cpp
