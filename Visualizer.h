//File to write Sols.

//for visualizing on Paraview
class ParaviewWriter {


  ParaviewWriter ParaviewWriter();

  void WriteOneTimeStep(const char* filename);

  void WriteAllTimeSteps(const char* filename);


  ~ParaviewWriter();
};
