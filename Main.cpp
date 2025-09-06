// One-way coupled analytical Fluid-Structure Riemann Problem (FSRP)
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include "Visualizer.h"
#include "RankineHugoniot.h"
using namespace std;


int main(){

  //! PARAMETERS

  //mesh
  double xmin = 0.0; //m
  double xmax = 5.0; //m
  double ptnum = 1.0e3;

  //initial conditions
  //fluid
  double rho_f = 1.3; //kg/m^3
  double vel_f = 0.0; //m/s
  double pressure_f = 1.0e5; //Pa
  double gamma = 1.4; //specific heat ratio

  //piston
  double rho_p = 0.0;
  double vel_p = 10.0;  // m/s
  double pressure_p = 0.0;

  //wall
  double rho_wall = 0.0;
  double vel_wall = 0.0;  // m/s
  double pressure_wall = 0.0;
  double wall_loc = 0.8*xmax; //wall starts at 80% of xmax

  //time(s)
  double t_final = 0.5;
  double t_samples = 1.0e3;

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
  double dt = t_final / t_samples;
  for (int n=0;n<=t_samples;n++)
    time.push_back(dt*n);
  

  //! DATA ALLOCATION
  //state variables
  vector<double> rho((int)xcoords.size());
  vector<double> vel((int)xcoords.size());
  vector<double> pressure((int)xcoords.size());
  //loc. of waves (piston + wall + shock)
  double shock_loc,piston_loc;
  vector<double> shock((int)time.size());
  vector<double> piston((int)time.size());


  //! COMPUTE INTERMEDIATE STATES -- Rankine-Hugoniot Jump Conditions
  //Note: assuming perfect gas
  //Note: from FSI Notes by Dr. Wang (found on Canvas)

  RH_BASE RH_Jump(gamma);
  array<double,3> shock_states = RH_Jump.ComputeNewStates(rho_f,vel_f,pressure_f,vel_p);
  double rho_star = shock_states[0]; 
  double pressure_star = shock_states[1]; 
  double vel_shock = shock_states[2];

  //! MAIN LOOP
  for (int t=0;t<=t_samples;t++){

    //compute new distance of piston + shock
    piston_loc = time[t] * vel_p;
    piston[t] = piston_loc;
    //piston[t] /= xmax;
    shock_loc = time[t] * vel_shock;
    shock[t] = shock_loc;
    //shock[t] /= xmax;
    
    for (int n=0;n<(int)xcoords.size();n++){ //iterating through all pts
      x = xcoords[n];
      
      if (x<=piston_loc){ //pts in piston
        rho[n] = rho_p; vel[n] = vel_p; pressure[n] = pressure_p;
      }
      else if (x>= wall_loc){ //pts in wall
        rho[n] = rho_wall; vel[n] = vel_wall;  pressure[n] = pressure_wall;
      }
      else { //pts in between wall + piston
       
        if (x<shock_loc){ //intermediate state
          rho[n] = rho_star; vel[n] = vel_p;  pressure[n] = pressure_star;
        }
        else if (x>shock_loc){ //undisturbed state
          rho[n] = rho_f; vel[n] = vel_f;  pressure[n] = pressure_f;
        }  
        else //error handling
          cout<<"Error: some pts left unassigned in between piston and wall!\n";
      }

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
  //Storing all time-steps of state variables in one file (.pvd)
  const char* file = "Results/Results.pvd";
  Visual.WriteAllTimeSteps(file,iter_visuals);

  //Writing all wave locations
  const char* file_contact = "Results/Contact_Locs.vtu";
  const char* file_shock = "Results/Shock_Locs.vtu";

  Visual.WriteWaveLoc(file_contact,piston,time);
  Visual.WriteWaveLoc(file_shock,shock,time);


  return 0;
}
