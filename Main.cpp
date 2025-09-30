// One-way coupled analytical Fluid-Structure Riemann Problem (FSRP)
#include <iostream>
#include <string> 
#include <vector>
#include <cmath>
#include <algorithm>
#include <array>

#include "Visualizer.h"
#include "RankineHugoniot.h"
using namespace std;


int main(){

  //! PARAMETERS

  //mesh
  double xmin = -2.0; //m
  double xmax = 5.0; //m
  double ptnum = 9.0e3;

  //initial conditions

  //fluid 1
  double rho_f1 = 1.3; //kg/m^3
  double vel_f1 = 0.0; //m/s
  double pressure_f1 = 1.0e5; //Pa
  double gamma1 = 1.4; //specific heat ratio

  //fluid 2
  double rho_f2 = 1.3; //kg/m^3
  double vel_f2 = 0.0; //m/s
  double pressure_f2 = 1.0e5; //Pa
  double gamma2 = 1.4; //specific heat ratio

  //piston
  double vel_p = 10.0;  // m/s
  double piston_loc = 0.0;

  //wall
  double rho_wall = 0.0;
  double vel_wall = 0.0;  // m/s
  double pressure_wall = 0.0;
  double wall_loc = xmax; 

  //time(s)
  double t_final = 0.005;
  //double t_final = 0.014;
  //double t_final = 0.033;
  //double t_final = 0.50;
  double t_samples = 1.0e2;

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
  //state variables (primitive vars.)
  vector<double> rho((int)xcoords.size());
  vector<double> vel((int)xcoords.size());
  vector<double> pressure((int)xcoords.size());
  double shock_loc,vel_s; //location + velocity of shock wave
  //location list of shock and piston
  vector<double> shock((int)time.size());
  vector<double> piston((int)time.size());
  array<double,3> left_state{rho_f2,vel_f2,pressure_f2}; //initial conditions are of fluid
  array<double,3> right_state{rho_f2,vel_f2,pressure_f2};
  double t_contact = 0.0;
  double piston_loc_wall_collision; //pos. piston and shock/wall contact each other


  //TODO: add limit to dt to prevent shock overshooting wall/piston by some max distance
  //TODO: add function to compute rarefaction fan

  RH_BASE RH_Jump(gamma2); //for handling shock wave
  array<double,3> shock_states;

  RFan_BASE RarFan(gamma1); //for handling rarefaction fan
  array<double,2> rarfan_loc; //index: [head,tail]
  array<double,3> rarfan_states; //state vars. in fan
  array<double,2> rarfan_stars; //state vars. affected by fan

  //NOTE: relative time vs. absolute time

  //! MAIN LOOP
  for (int t=0;t<=t_samples;t++){

    //! COMPUTE INTERMEDIATE STATES 

    // SHOCK - RANKINE-HUGONIOT JUMP CONDITIONS
    //1st time -- compute intermediate state
    if (t==0) {
      shock_states = RH_Jump.ComputeShockStates(rho_f2,vel_f2,pressure_f2,vel_p);
      vel_s = shock_states[2];
      left_state[0] = shock_states[0];left_state[1] = vel_p; left_state[2] = shock_states[1]; //new state of fluid affected by shock
    }

    //compute new distance of piston + shock
    piston_loc = time[t] * vel_p;
    //note: positive vel_s refers to shock moving towards away from the piston and towards the wall
    shock_loc = (vel_s>0) ? piston_loc_wall_collision + (time[t]-t_contact) * vel_s : wall_loc + (time[t]-t_contact) * vel_s;

    if (shock_loc >= wall_loc){ //SHOCK-WALL CASE
      shock_states = RH_Jump.ComputeShockStates(left_state[0],left_state[1],left_state[2],vel_wall);
      right_state[0] = shock_states[0];right_state[1] = vel_wall; right_state[2] = shock_states[1]; //new state of fluid affected by shock

      //correcting shock postion 
      t_contact = time[t] - (shock_loc-wall_loc)/(vel_s); //time when shock contacts wall
      piston_loc_wall_collision = t_contact * vel_p; //pos. piston is when shock contacted wall -- used for following if statement

      shock_loc = wall_loc + shock_states[2]*(time[t]-t_contact); //corrected position
      vel_s = shock_states[2]; //updating shock velocity

    }
    if ((shock_loc <= piston_loc) && (t!=0)){ //SHOCK-PISTON CASE
      shock_states = RH_Jump.ComputeShockStates(right_state[0],right_state[1],right_state[2],vel_p);
      left_state[0] = shock_states[0];left_state[1] = vel_p; left_state[2] = shock_states[1]; //new state of fluid affected by shock

      //correcting shock postion
      //t_contact += (wall_loc - piston_loc_wall_collision) / (vel_p-vel_s);
      t_contact = (shock_loc - piston_loc - vel_s*time[t] + vel_p*time[t]) / (vel_p-vel_s);
      piston_loc_wall_collision = t_contact * vel_p; //pos. piston is when shock contacted wall -- used for following if statement
      shock_loc = piston_loc + shock_states[2]*(time[t]-t_contact); //corrected position
      vel_s = shock_states[2]; //updating shock velocity
    }

    shock[t] = shock_loc;
    piston[t] = piston_loc;

    //RAREFACTION FAN LOCATION
    rarfan_loc = RarFan.ComputeRareFanLocs(rho_f1,vel_f1,pressure_f1,vel_p,time[t]);

    //! ASSIGNING STATE VARS. TO FIELD
    for (int n=0;n<(int)xcoords.size();n++){ //iterating through all pts
      x = xcoords[n];

      //pts in fluid 1 (rarefaction)
      if (x<=piston_loc){ 
        if (x>=rarfan_loc[0] && x<=rarfan_loc[1]){ //pts in rarefaction fan
          rarfan_states = RarFan.ComputeRareFanStates(x,time[t],rho_f1,vel_f1,pressure_f1);
          rho[n] = rarfan_states[0]; vel[n] = rarfan_states[1]; pressure[n] = rarfan_states[2];
        }
        else if (x>rarfan_loc[1]){ //pts in between piston and rarefaction fan
          rarfan_stars = RarFan.ComputeStarVals(rho_f1,vel_f1,pressure_f1,vel_p);
          rho[n] = rarfan_stars[0]; vel[n] = vel_p; pressure[n] = rarfan_stars[1];
        }
        else{  //undisturbed fluid
          rho[n] = rho_f1; vel[n] = vel_f1; pressure[n] = pressure_f1;
        }
      }

      //pts in wall
      else if (x>= wall_loc){
        rho[n] = rho_wall; vel[n] = vel_wall;  pressure[n] = pressure_wall;
      }

      //pts in between wall + piston
      else { 
       
        if (x<shock_loc){ //left state of shock
          rho[n] = left_state[0]; vel[n] = left_state[1];  pressure[n] = left_state[2];
        }
        else if (x>shock_loc){ //right state of shock
          rho[n] = right_state[0]; vel[n] = right_state[1];  pressure[n] = right_state[2];
        }  
        else //error handling
          cout<<"Error: some pts left unassigned in between piston and wall!\n";
      }

    //print out results and save for visualization
    }
    if (t%iterout == 0){
      time_step = to_string(time[t]);
      time_step += "s";
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
