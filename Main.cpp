// One-way coupled analytical Fluid-Structure Riemann Problem (FSRP)
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include "Visualizer.h"
using namespace std;


int main(){

  //! PARAMETERS

  //mesh
  double xmin = 0.0; //m
  double xmax = 1.0; //m
  double ptnum = 10.0;

  //initial conditions
  //fluid
  double rho_f = 1.3; //kg/m^3
  double vel_f = 0.0; //m/s
  double pressure_f = 1.0e5; //Pa
  double gamma = 1.4; //specific heat ratio

  //piston
  double rho_p = 0.0;
  double vel_p = 5.0;  // m/s
  double pressure_p = 0.0;

  //time(s)
  double t_max = 1.0;
  double t_samples = 10.0;

  //visualization
  ParaviewWriter Visual;
  int iterout = 2;
  vector<string> iter_visuals; //stores all iter. out #s for writing .pvd file
  string name,time_step;

  //! MESH + TIME GENERATION
  //xcoords
  vector<double> xcoords;
  double dx = (xmax-xmin) / ptnum;
  double x; //current xcoord
  for (int n=0;n<=ptnum;n++)
    xcoords.push_back(xmin + dx*n);
  //time
  vector<double> time;
  double dt = t_max / t_samples;
  for (int n=0;n<=t_samples;n++)
    time.push_back(dt*n);
  

  //! DATA ALLOCATION
  //state variables
  vector<double> rho((int)xcoords.size());
  vector<double> vel((int)xcoords.size());
  vector<double> pressure((int)xcoords.size());
  //loc. of waves (piston + shock)
  double shock_loc,piston_loc;

  //! COMPUTE INTERMEDIATE STATES -- Rankine-Hugoniot Jump Conditions
  //Note: assuming perfect gas
  double alpha = 2.0 / (gamma+1.0)*rho_f;
  double beta = pressure_f * ( (gamma-1.0) / (gamma+1.0) );
  double phi = alpha / pow(vel_p-vel_f,2.0);
  
  double pressure_star = ( 1.0 + sqrt(4.0*phi*(pressure_f+beta)+1.0) ) / (2.0*phi);
  pressure_star += pressure_f; //pressure
  double rho_star = (pressure_star/pressure_f) + ((gamma-1.0)/(gamma+1.0)) /
  ( (pressure_star/pressure_f) * ((gamma-1.0)/(gamma+1.0)) );
  rho_star *= rho_f; //density

  double vel_shock = (rho_f*vel_f - rho_star*vel_p) / (rho_f-rho_star); //shock speed

  //! MAIN LOOP
  for (int t=0;t<t_samples;t++){

    //compute new distance of piston + shock
    piston_loc = time[t] * vel_p;
    shock_loc = time[t] * vel_shock;
    
    for (int n=0;n<(int)xcoords.size();n++){ //iterating through all pts
      x = xcoords[n];
      // pts in piston
      if (x<=piston_loc){
        rho[n] = rho_p; vel[n] = vel_p; pressure[n] = pressure_p;
      }
      //pts in intermediate region 
      else if (x>piston_loc && x<=shock_loc){
        rho[n] = rho_star; vel[n] = vel_p;  pressure[n] = pressure_star;
      }
      else if (x>shock_loc){
        rho[n] = rho_f; vel[n] = vel_f;  pressure[n] = pressure_f;
      } 
      else //error handling
        cout<<"Error!\n";

    //print out results and save for visualization
    }
    if (t%iterout == 0){
      time_step = to_string(t);
      iter_visuals.push_back(time_step);
      name = "Results/Time";
      name += time_step;
      name += ".vtu";
      const char* filename = name.c_str();
      Visual.WriteOneTimeStep(filename,xcoords,rho,vel,pressure);
    }
  }

  //! POST-PROCESSING
  //Storing all time-steps in one file (.pvd)
  const char* file = "Results/Results.pvd";
  Visual.WriteAllTimeSteps(file,iter_visuals);


  return 0;
}
