//File to write Sols.
#ifndef VISUALIZER_H
#define VISUALIZER_H

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>

using namespace std;

//for visualizing on Paraview
class ParaviewWriter {

  public:
  ParaviewWriter();

  void WriteOneTimeStep(const char* &filename,vector<double> &xcoords,vector<double> &density,vector<double> &velocity,vector<double> &pressure);

  void WriteAllTimeSteps(const char* &filename,vector<string> &iter_visuals);

  void WriteWaveLoc(const char* &filename,vector<double> &wave_pos,vector<double> &times);

  ~ParaviewWriter();
};
#endif 
