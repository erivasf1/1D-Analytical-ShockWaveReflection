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
  double xmin = -4.0; //m
  double xmax = 5.0; //m
  double ptnum = 1.0e3;
  //double ptnum = 10;

  //scenario
  //note: 1: overlapped region assigned to fluid 1 (only shock in fluid 2)
  //      2: overlapped region left unassigned/new fluid (shock in fluid 1&2)
  int scenario = 1;
  bool top_change{false};

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

  //fluid 3
  double rho_f3 = 1.3; //kg/m^3
  double vel_f3 = 0.0; //m/s
  double pressure_f3 = 1.0e5; //Pa
  double gamma3 = 1.4; //specific heat ratio

  //piston
  double vel_p = 10.0;  // m/s
  double piston_loc = 0.0;

  //wall
  double vel_wall = 0.0;  // m/s
  double wall_loc = 2.0; 

  //time(s)
  //double t_final = 0.005;
  //double t_final = 0.025;
  double t_final = 0.03;
  //double t_final = 0.50;
  //double t_samples = 1.0e3;
  double t_samples = 5.0;
  double t_topchange = (wall_loc-piston_loc) / vel_p; //time when piston interpenetrates wall

  //visualization
  ParaviewWriter Visual;
  int iterout = 1;
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

  //shock + rarefacetion fan vars.
  //vector<double> shock((int)time.size());
  //vector<double> piston((int)time.size());
  [[maybe_unused]] array<double,3> left_state1{rho_f1,vel_f1,pressure_f1}; //initial conditions are of fluid3
  [[maybe_unused]] array<double,3> right_state1{rho_f1,vel_f1,pressure_f1};
  array<double,3> left_state2{rho_f2,vel_f2,pressure_f2}; //initial conditions are of fluid3
  array<double,3> right_state2{rho_f2,vel_f2,pressure_f2};
  array<double,3> left_state3{rho_f3,vel_f3,pressure_f3}; //initial conditions are of fluid3
  array<double,3> right_state3{rho_f3,vel_f3,pressure_f3};

  double shock_loc1;//vel_s1; //location + velocity of shock wave
  double shock_loc2;//vel_s2; //location + velocity of shock wave
  double shock_loc3;//vel_s3; //location + velocity of shock wave

  array<double,3> shock_state1,shock_state3;

  //TODO: IMPROVE ANALYTICAL SOL. TO HANDLE BOTH SCENARIOS INSTEAD OF SIMPLE FSRP!!!

  RH_BASE RH_Jump1(gamma1,rho_f1,vel_f1,pressure_f1,vel_p,vel_wall,piston_loc,wall_loc); //for handling shock wave
  RH_BASE RH_Jump2(gamma2,rho_f2,vel_f2,pressure_f2,vel_p,vel_wall,piston_loc,wall_loc); //for handling shock wave
  RH_BASE RH_Jump3(gamma3,rho_f3,vel_f3,pressure_f3,vel_p,vel_wall,piston_loc,wall_loc); //for handling shock wave

  RFan_BASE RarFan(gamma1); //for handling rarefaction fan -- only for fluid 1
  array<double,2> rarfan_loc; //index: [head,tail]
  array<double,3> rarfan_states; //state vars. in fan
  array<double,2> rarfan_stars; //state vars. affected by fan

  //! MAIN LOOP
  for (int t=0;t<=t_samples;t++){

    piston_loc = time[t] * vel_p; //piston loc.

    top_change = (time[t] >= t_topchange) ? true : false;

    //FLUID 2 (piston-dependent) & 3 (piston-dependent) SHOCKS
    if (top_change == false) //only computing shock in fluid 2 if piston didn't interpenetrate wall
      RH_Jump2.ComputeShockResults(time[t],left_state2,right_state2,shock_loc2);
    else{
      shock_state3 = RH_Jump3.ComputeShockStates(right_state3[0],right_state3[1],right_state3[2],vel_p);
      shock_loc3 = ( shock_state3[2] * (time[t]-t_topchange)) + wall_loc;
      left_state3[0] = shock_state3[0]; left_state3[1] = vel_p; left_state3[2] = shock_state3[1];
    }


    //FLUID 1 RAREFACTION FAN & SHOCK (scenario&piston-dependent)
    rarfan_loc = RarFan.ComputeRareFanLocs(rho_f1,vel_f1,pressure_f1,vel_p,time[t]);
    rarfan_stars = RarFan.ComputeStarVals(rho_f1,vel_f1,pressure_f1,vel_p); //affected fluid in between rarefaction fan & piston
    if ( (top_change == true) && (scenario==2) ){ //shock in fluid 1 case 
      left_state1[0]=rarfan_stars[0];left_state1[1]=vel_p,left_state1[2]=rarfan_stars[1];
      shock_state1 = RH_Jump1.ComputeShockStates(left_state1[0],left_state1[1],left_state1[2],vel_wall);
      shock_loc1 = wall_loc + shock_state1[2] * (time[t] - t_topchange);
      right_state1[0] = shock_state1[0]; right_state1[1] = vel_wall; right_state1[2] = shock_state1[1];
    }


    //! ASSIGNING STATE VARS. TO FIELD
    for (int n=0;n<(int)xcoords.size();n++){ //iterating through all pts
      x = xcoords[n];

      //pts in FLUID 1 
      if (x<=piston_loc){ 
        if ( (top_change == true) && (scenario==2) ){ //scenario 2 after piston interpenetration
          if (x<shock_loc1){
            rho[n] = left_state1[0]; vel[n] = left_state1[1]; pressure[n] = left_state1[2];
          }      
          else if ((x>shock_loc1) && (x<wall_loc)) {
            rho[n] = right_state1[0]; vel[n] = right_state1[1]; pressure[n] = right_state1[2];
          }
          else if ((x>wall_loc) && (x<piston_loc) ){ //inactive region (setting to initial condition)
            rho[n] = rho_f1; vel[n] = vel_f1; pressure[n] = pressure_f1;
            //rho[n] = 0.0; vel[n] = 0.0; pressure[n] = 0.0;
          }
        }
        else { //before piston interpenetration or scenario 1
          if (x>=rarfan_loc[0] && x<=rarfan_loc[1]){ //pts in rarefaction fan
            rarfan_states = RarFan.ComputeRareFanStates(x,time[t],rho_f1,vel_f1,pressure_f1);
            rho[n] = rarfan_states[0]; vel[n] = rarfan_states[1]; pressure[n] = rarfan_states[2];
          }
          else if (x>rarfan_loc[1]){ //pts in between piston and rarefaction fan
            rho[n] = rarfan_stars[0]; vel[n] = vel_p; pressure[n] = rarfan_stars[1];
          }
          else{  //undisturbed fluid
            rho[n] = rho_f1; vel[n] = vel_f1; pressure[n] = pressure_f1;
          }

        }
      }

      //pts in fluid 2 (fluid 2 can only exist when piston doesn't interpenetrate wall)
      else if ( (top_change == false) && (x < wall_loc) ) {
        if (x<shock_loc2){
          rho[n] = left_state2[0]; vel[n] = left_state2[1]; pressure[n] = left_state2[2];
        }
        else if (x>shock_loc2){
          rho[n] = right_state2[0]; vel[n] = right_state2[1]; pressure[n] = right_state2[2];
        }
        else
          cerr<<"Error: Pts left unassigned in fluid 2!"<<endl;
      }

      //pts in fluid 3
      else if ((x>wall_loc) && (x>piston_loc)){
        if (top_change == true){ //piston interpenetration scenario -- causing shock
          if (x<shock_loc3){
            rho[n] = left_state3[0]; vel[n] = left_state3[1]; pressure[n] = left_state3[2];
          }
          else {
            rho[n] = right_state3[0]; vel[n] = right_state3[1]; pressure[n] = right_state3[2];
          }
        }
        else {
          rho[n] = left_state3[0]; vel[n] = left_state3[1]; pressure[n] = left_state3[2];
        }
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
  //const char* file_contact = "Results/Contact_Locs.vtu";
  //const char* file_shock = "Results/Shock_Locs.vtu";

  //Visual.WriteWaveLoc(file_contact,piston,time);
  //Visual.WriteWaveLoc(file_shock,shock,time);


  return 0;
}
