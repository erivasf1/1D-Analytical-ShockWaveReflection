//Class computing intermediate states via Rankine-Hugoniot Jump Conditions
#ifndef RANKINEHUGONIOT_H
#define RANKINEHUGONIOT_H

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace std;

class RH_BASE {
 
  //NOTE: using perfect gas relations
  double gamma;

  public:
  RH_BASE(double g);

  double ComputeAlpha(double rho_f);
  double ComputeBeta(double pressure_f);
  double ComputePhi(double vel_f,double vel_p,double alpha);

  array<double,3> ComputeNewStates(double rho_f,double vel_f,double pressure_f,double vel_p);

  ~RH_BASE();

};
#endif 
