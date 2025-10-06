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
  double gamma; //specific heat-ratio
  double rho_f,vel_f,pressure_f; //initial fluid states
  double vel_p,vel_w; //piston vel. & wall vel.
  double piston_loc,wall_loc; //piston & wall initial states

  public:
  RH_BASE(double g,double r,double v,double press,double vp,double vw,double ploc,double wloc);

  double ComputeAlpha(double rho);
  double ComputeBeta(double pressure);
  double ComputePhi(double vel,double alpha,double vel_wall);

  array<double,3> ComputeShockStates(double rho,double vel,double pressure,double vel_wall); //inputs=undisturbed fluid states

  void ComputeShockResults(double t_target,array<double,3> &left_state,array<double,3> &right_state,double &shock_loc); //assigns vals. for the left and right state, & the shock loc.

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
