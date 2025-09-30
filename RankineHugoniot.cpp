//Defs for Rankine-Hugoniot Class
#include "RankineHugoniot.h"
//----------------- RH_BASE CLASS DEFS. ------------------------------
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
array<double,3> RH_BASE::ComputeShockStates(double rho_f,double vel_f,double pressure_f,double vel_p){

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
void RH_BASE::AssignIntermediateStates(vector<double> &rho,vector<double> &vel,vector<double> &pressure,vector<double> &xcoords,double shock_loc,array<double,3> &left_state,array<double,3> &right_state,bool dir,double piston_loc,double wall_loc){

  double x;
  //NOTE: dir refers to the dir. of shock
  for (int n=0;n<(int)xcoords.size();n++){ //iterating through all pts
    x = xcoords[n];

    if ((x<=piston_loc) || (x>= wall_loc)) //pts in piston or in wall case
      continue;

    //+x dir. of shock case
    if (dir == true){
      if (x<shock_loc){
        rho[n] = left_state[0]; vel[n] = left_state[1];  pressure[n] = left_state[2];
      }
      else{
        rho[n] = right_state[0]; vel[n] = right_state[1];  pressure[n] = right_state[2];
      }
    }

    //-x dir. of shock case
    else {
      if (x<shock_loc){
        rho[n] = right_state[0]; vel[n] = right_state[1];  pressure[n] = right_state[2];
      }
      else{
        rho[n] = left_state[0]; vel[n] = left_state[1];  pressure[n] = left_state[2];
      }

    }

  }

  return;
}
//--------------------------------------------------------------------
RH_BASE::~RH_BASE(){}
//--------------------------------------------------------------------

//----------------- RFan_BASE CLASS DEFS. ------------------------------

//--------------------------------------------------------------------
RFan_BASE::RFan_BASE(double g) : gamma(g) {}
//--------------------------------------------------------------------
array<double,2> RFan_BASE::ComputeRareFanLocs(double rho_f,double vel_f,double pressure_f,double vel_p,double t){

  array<double,2> starvals = ComputeStarVals(rho_f,vel_f,pressure_f,vel_p);

  double a = ComputeSoundSpeed(rho_f,pressure_f);
  double a_star = ComputeSoundSpeed(starvals[0],starvals[1]);

  double left_end = (vel_f - a) * t; double right_end = (vel_p - a_star) * t;
  array<double,2> locs{left_end,right_end};

  return locs;

}
//--------------------------------------------------------------------
array<double,2> RFan_BASE::ComputeStarVals(double rho_f,double vel_f,double pressure_f,double vel_p){

  //intermediate pressure
  double a = ComputeSoundSpeed(rho_f,pressure_f);
  double p_star = -((gamma-1.0)/(2.0*a)) * (vel_p-vel_f) + 1.0;
  p_star = pressure_f * pow(p_star,(2.0*gamma)/(gamma-1.0));

  //intermediate density
  double rho_star = rho_f* pow(p_star/pressure_f,1/gamma);

  array<double,2> starvals{rho_star,p_star};
  return starvals;

}
//--------------------------------------------------------------------
array<double,3> RFan_BASE::ComputeRareFanStates(double xcoord,double t,double rho_f,double vel_f,double pressure_f){

  double alpha = xcoord / t;
  double a = ComputeSoundSpeed(rho_f,pressure_f);

  //velocity
  double v = ( (gamma-1.0)*vel_f + 2.0*(a+alpha) ) / (gamma+1.0);
  //density
  double rho = ( pow(rho_f,gamma) * pow(v-alpha,2.0) ) / (gamma*pressure_f);
  rho = pow(rho,1.0/(gamma-1.0));
  //pressure
  double p = pressure_f * pow(rho/rho_f,gamma);

  array<double,3> states{rho,v,p};

  return states;

}
//--------------------------------------------------------------------
double RFan_BASE::ComputeSoundSpeed(double rho,double p){

  double res = sqrt( (gamma*p) / rho );
  return res;

}
//--------------------------------------------------------------------
RFan_BASE::~RFan_BASE(){}
//--------------------------------------------------------------------



