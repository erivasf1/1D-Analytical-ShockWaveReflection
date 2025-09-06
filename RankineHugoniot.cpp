//Defs for Rankine-Hugoniot Class
#include "RankineHugoniot.h"
//--------------------------------------------------------------------
RH_BASE::RH_BASE(double g) : gamma(g) {}
//--------------------------------------------------------------------
double RH_BASE::ComputeAlpha(double rho_f){

  double res = 2.0 / ( (gamma+1.0)*rho_f );
  return res;

}
//--------------------------------------------------------------------
double RH_BASE::ComputeBeta(double pressure_f){

  double res = pressure_f * ( (gamma-1.0) / (gamma+1.0) );
  return res;

}
//--------------------------------------------------------------------
double RH_BASE::ComputePhi(double vel_f,double vel_p,double alpha){

  double res = alpha / pow(vel_p-vel_f,2.0);
  return res;

}
//--------------------------------------------------------------------
array<double,3> RH_BASE::ComputeNewStates(double rho_f,double vel_f,double pressure_f,double vel_p){

  double alpha = ComputeAlpha(rho_f);
  double beta = ComputeBeta(pressure_f);
  double phi = ComputePhi(vel_f,vel_p,alpha);

  double pressure_star = ( 1.0 + sqrt(4.0*phi*(pressure_f+beta)+1.0) ) / (2.0*phi);
  pressure_star += pressure_f; //pressure
  double rho_star = ( (pressure_star/pressure_f) + ((gamma-1.0)/(gamma+1.0)) ) /
  ( (pressure_star/pressure_f) * ((gamma-1.0)/(gamma+1.0)) + 1.0 );
  rho_star *= rho_f; //density

  double vel_shock = (rho_f*vel_f - rho_star*vel_p) / (rho_f-rho_star); //shock speed
  array<double,3> sols{rho_star,pressure_star,vel_shock};

  return sols;
}
//--------------------------------------------------------------------
RH_BASE::~RH_BASE(){}
//--------------------------------------------------------------------
