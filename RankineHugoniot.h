//Class computing intermediate states via Rankine-Hugoniot Jump Conditions for shock and isentropic conditions for rarefaction fan
#ifndef RANKINEHUGONIOT_H
#define RANKINEHUGONIOT_H

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <array>

using namespace std;

class RH_BASE {
 
  //NOTE: using perfect gas relations
  double gamma;

  public:
  RH_BASE(double g);

  double ComputeAlpha(double rho_f);
  double ComputeBeta(double pressure_f);
  double ComputePhi(double vel_f,double vel_p,double alpha);

  array<double,3> ComputeShockStates(double rho_f,double vel_f,double pressure_f,double vel_p);

  void AssignIntermediateStates(vector<double> &rho,vector<double> &vel,vector<double> &pressure,vector<double> &xcoords,double shock_loc,array<double,3> &left_state,array<double,3> &right_state,bool dir,double piston_loc,double wall_loc);

  ~RH_BASE();

};

class RFan_BASE {

  double gamma;

  //assuming isentropic conditions & perfect gas relations
  public:
  RFan_BASE(double g);

  array<double,2> ComputeRareFanLocs(double rho_f,double vel_f,double pressure_f,double vel_p,double t); //returns the left(head) and right(tail) end of the rarefaction fan;index: (left,right)
  array<double,2> ComputeStarVals(double rho_f,double vel_f,double pressure_f,double vel_p);

  array<double,3> ComputeRareFanStates(double xcoord,double t,double rho_f,double vel_f,double pressure_f); //[0]:rho [1]:vel [2]:pressure

  double ComputeSoundSpeed(double rho,double p);


  ~RFan_BASE();

};

#endif 
