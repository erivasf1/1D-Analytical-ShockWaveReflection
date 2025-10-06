//Defs for Rankine-Hugoniot Class
#include "RankineHugoniot.h"
//----------------- RH_BASE CLASS DEFS. ------------------------------
//--------------------------------------------------------------------
RH_BASE::RH_BASE(double g,double r,double v,double press,double vp,double vw,double ploc,double wloc) : gamma(g), rho_f(r), vel_f(v), pressure_f(press), vel_p(vp), vel_w(vw), piston_loc(ploc), wall_loc(wloc) {}
//--------------------------------------------------------------------
double RH_BASE::ComputeAlpha(double rho){

  double res = 2.0 / ( (gamma+1.0)*rho );
  return res;

}
//--------------------------------------------------------------------
double RH_BASE::ComputeBeta(double pressure){

  double res = pressure * ( (gamma-1.0) / (gamma+1.0) );
  return res;

}
//--------------------------------------------------------------------
double RH_BASE::ComputePhi(double vel,double alpha,double vel_wall){

  double res = alpha / pow(vel_wall-vel,2.0);
  return res;

}
//--------------------------------------------------------------------
array<double,3> RH_BASE::ComputeShockStates(double rho,double vel,double pressure,double vel_wall){
  //note: vel_wall is either the vel. of the moving piston or stationary wall

  double alpha = ComputeAlpha(rho);
  double beta = ComputeBeta(pressure);
  double phi = ComputePhi(vel,alpha,vel_wall);

  double pressure_star = ( 1.0 + sqrt(4.0*phi*(pressure+beta)+1.0) ) / (2.0*phi);
  pressure_star += pressure; //pressure

  double rho_star = ( (pressure_star/pressure) + ((gamma-1.0)/(gamma+1.0)) ) /
  ( (pressure_star/pressure) * ((gamma-1.0)/(gamma+1.0)) + 1.0 );
  rho_star *= rho; //density

  double vel_shock = (rho*vel - rho_star*vel_wall) / (rho-rho_star); //shock speed

  array<double,3> sols{rho_star,pressure_star,vel_shock}; //computed shock states

  return sols;
}
//--------------------------------------------------------------------
void RH_BASE::ComputeShockResults(double t_target,array<double,3> &left_state,array<double,3> &right_state,double &shock_loc){

  if (t_target <= 0.0)
    return;

  double t_contactw = 0.0; double t_contactp = 0.0; //collision times of shock with wall & piston,respectively
  array<double,3> shock_states; //{rho_star,pressure_star,vel_shock}

  left_state[0] = rho_f; left_state[1] = vel_f; left_state[2] = pressure_f; //!< resetting primitive vars. from previous time computation
  right_state[0] = rho_f; right_state[1] = vel_f; right_state[2] = pressure_f;

  for (int n=0;n<1e4;n++){
    //Step 1: compute new left state from RH_JUMP conds.
    //shock_states = ComputeShockStates(rho_f,vel_f,pressure_f,vel_p);
    shock_states = ComputeShockStates(right_state[0],right_state[1],right_state[2],vel_p);
    //left_state(shock_states[0],vel_p,shock_states[1]);
    left_state[0] = shock_states[0]; left_state[1] = vel_p; left_state[2] = shock_states[1];

    //Step 2: compare if target time is greater than time of shock-wall collision
    // if not, compute position of shock and piston and break. if yes, go to step 3
    t_contactw = ( (wall_loc - piston_loc) / shock_states[2] ) + t_contactp;
    if (t_contactw >= t_target){
      shock_loc = piston_loc + shock_states[2] * (t_target-t_contactp);
      break;
    }
    piston_loc = vel_p * t_contactw; //upating piston loc.
    
    //Step 3: compute new right state from RH_JUMP conds.
    shock_states = ComputeShockStates(left_state[0],left_state[1],left_state[2],vel_w);
    //right_state{shock_states[0],vel_w,shock_states[1]};
    right_state[0] = shock_states[0]; right_state[1] = vel_w; right_state[2] = shock_states[1];

    //Step 4: compare if target time is greater than time of piston collision
    // if not, compute position of shock and piston and break. if yes, go to step 2
    t_contactp = (wall_loc - piston_loc - shock_states[2]*t_contactw) / (vel_p-shock_states[2]); //TODO: need to look into this further!!!
    if (t_contactp >= t_target){
      shock_loc = wall_loc + shock_states[2] * (t_target-t_contactw);
      break;
    }
    piston_loc = vel_p * t_contactp; //updating piston loc.
  
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



