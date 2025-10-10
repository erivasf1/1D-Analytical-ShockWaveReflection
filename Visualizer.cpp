//Visualizer Defs.
#include "Visualizer.h"

//Paraview Class Defs.
//-------------------------------------------------
ParaviewWriter::ParaviewWriter(){}
//-------------------------------------------------
void ParaviewWriter::WriteOneTimeStep(const char* &filename,vector<double> &xcoords,vector<double> &density,vector<double> &velocity,vector<double> &pressure){
  //Note: follows .vtu format (unstructured)

  ofstream myfile(filename);

  if (!myfile){ //checking if file opened successfully
    cerr<<"Error: Could Not Open File "<<filename<<endl;
    return;
  }

  //TITLE
  myfile<<"<?xml version=\"1.0\"?>"<<endl;
  myfile<<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">"<<endl;
  myfile<<"  <UnstructuredGrid>"<<endl;
  myfile<<"    <Piece NumberOfPoints=\""<<(int)xcoords.size()<<"\" NumberOfCells=\""<<(int)xcoords.size()-1<<"\">"<<endl;
  myfile<<endl;
  
  //COORDS. INFO.
  myfile<<"      <!-- Point coordinates -->"<<endl;
  myfile<<"      <Points>"<<endl;
  myfile<<"        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">"<<endl;
  for (int n=0;n<(int)xcoords.size();n++)
    myfile<<"        "<<xcoords[n]<<" 0.0"<<" 0.0"<<endl;    

  myfile<<"        </DataArray>"<<endl;
  myfile<<"      </Points>"<<endl;
  myfile<<endl;

  int count = 0; //for ensuring 6 columns of data per data type (e.g. coords)
  //CELLS (CONNECTIVITIES) INFO.
  myfile<<"      <!-- Connectivity (lines between points) -->"<<endl;
  myfile<<"      <Cells>"<<endl;
  myfile<<"        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">"<<endl;
  for (int n=0;n<(int)xcoords.size()-1;n++)
    myfile<<"        "<<n<<" "<<n+1<<endl;

  myfile<<"        </DataArray>"<<endl;
  myfile<<"        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">"<<endl;
  myfile<<"        ";
  for (int n=1;n<(int)xcoords.size();n++){
    count++;
    myfile<<2*n<<" ";    
    if (count % 6 == 0){
      myfile<<endl;
      myfile<<"        ";
    }
  }
  count = 0;
  myfile<<endl;
  myfile<<"        </DataArray>"<<endl;
  myfile<<"        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">"<<endl;
  myfile<<"        ";
  for (int n=0;n<(int)xcoords.size();n++){
    count++;
    myfile<<3<<" ";
    if (count % 6 == 0){
      myfile<<endl;
      myfile<<"        ";
    }
  }
  count = 0;
  myfile<<endl;
  myfile<<"        </DataArray>"<<endl;
  myfile<<"      </Cells>"<<endl;

  //STATE VARIABLES INFO.
  myfile<<"      <!-- Fields at points -->"<<endl;
  myfile<<"      <PointData Scalars=\"scalars\">"<<endl;
  myfile<<"        <DataArray type=\"Float32\" Name=\"density\" format=\"ascii\">"<<endl;
  myfile<<"        ";
  for (int n=0;n<(int)density.size();n++){
    count++;
    myfile<<density[n]<<" ";
    if (count % 6 == 0){
      myfile<<endl;
      myfile<<"        ";
    }
  }
  count = 0;
  myfile<<endl;
  myfile<<"        </DataArray>"<<endl;
  myfile<<"        <DataArray type=\"Float32\" Name=\"velocity\" format=\"ascii\">"<<endl;
  myfile<<"        ";
  for (int n=0;n<(int)velocity.size();n++){
    count++;
    myfile<<velocity[n]<<" ";
    if (count % 6 == 0){
      myfile<<endl;
      myfile<<"        ";
    }
  }
  count = 0;
  myfile<<endl;
  myfile<<"        </DataArray>"<<endl;
  myfile<<"        <DataArray type=\"Float32\" Name=\"pressure\" format=\"ascii\">"<<endl;
  myfile<<"        ";
  for (int n=0;n<(int)pressure.size();n++){
    count++;
    myfile<<pressure[n]<<" ";
    if (count % 6 == 0){
      myfile<<endl;
      myfile<<"        ";
    }
  }
  count = 0;
  myfile<<endl;
  myfile<<"        </DataArray>"<<endl;
  myfile<<"      </PointData>"<<endl;

  //END OF FILE
  myfile<<"    </Piece>"<<endl;
  myfile<<"  </UnstructuredGrid>"<<endl;
  myfile<<"</VTKFile>"<<endl;

/* --Template
<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">
  <UnstructuredGrid>
    <Piece NumberOfPoints="5" NumberOfCells="4">

      <!-- Point coordinates -->
      <Points>
        <DataArray type="Float32" NumberOfComponents="3" format="ascii">
          0.0 0.0 0.0
          1.0 0.0 0.0
          2.0 0.0 0.0
          3.0 0.0 0.0
          4.0 0.0 0.0
        </DataArray>
      </Points>

      <!-- Connectivity (lines between points) -->
      <Cells>
        <DataArray type="Int32" Name="connectivity" format="ascii">
          0 1
          1 2
          2 3
          3 4
        </DataArray>
        <DataArray type="Int32" Name="offsets" format="ascii"> -> apply count here
          2 4 6 8
        </DataArray>
        <DataArray type="UInt8" Name="types" format="ascii"> -> apply count here
          3 3 3 3
        </DataArray>
      </Cells>

      <!-- Fields at points -->
      <PointData Scalars="scalars">
        <DataArray type="Float32" Name="density" format="ascii"> -> apply count here
          1.0 0.9 0.85 0.95 1.0
        </DataArray>
        <DataArray type="Float32" Name="velocity" format="ascii"> -> apply count here
          0.0 0.1 0.2 0.1 0.0
        </DataArray>
        <DataArray type="Float32" Name="pressure" format="ascii"> -> apply count here
          1.0 0.95 0.92 0.97 1.0
        </DataArray>
      </PointData>

    </Piece>
  </UnstructuredGrid>
</VTKFile>
  */
  return;
}
//-------------------------------------------------
void ParaviewWriter::WriteAllTimeSteps(const char* &filename,vector<string> &iter_visuals){

  ofstream myfile(filename);

  if (!myfile){ //checking if file opened successfully
    cerr<<"Error: Could Not Open File "<<filename<<endl;
    return;
  }

  //TITLE
  myfile<<"<?xml version=\"1.0\"?>"<<endl;
  myfile<<"<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">"<<endl;
  myfile<<"  <Collection>"<<endl;
  //Writing file names
  for (int n=0;n<(int)iter_visuals.size();n++)
    myfile<<"    <DataSet timestep=\""<<iter_visuals[n]<<"\""<<" part=\"0\""<<" file=\"Time"<<iter_visuals[n]<<".vtu\"/>"<<endl;

  myfile<<"  </Collection>"<<endl;
  myfile<<"</VTKFile>"<<endl;
/*
<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">
  <Collection>
    <DataSet timestep="0" part="0" file="solution_0.vtu"/>
    <DataSet timestep="1" part="0" file="solution_1.vtu"/>
    <DataSet timestep="2" part="0" file="solution_2.vtu"/>
    <DataSet timestep="3" part="0" file="solution_3.vtu"/>
    <!-- Add more timesteps here -->
  </Collection>
</VTKFile>
*/

  return;
}
//-------------------------------------------------
void ParaviewWriter::WriteWaveLoc(const char* &filename,vector<double> &wave_pos,vector<double> &times){

  ofstream myfile(filename);

  if (!myfile){ //checking if file opened successfully
    cerr<<"Error: Could Not Open File "<<filename<<endl;
    return;
  }

  //TITLE
  myfile<<"<?xml version=\"1.0\"?>"<<endl;
  myfile<<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">"<<endl;
  myfile<<"  <UnstructuredGrid>"<<endl;
  myfile<<"    <Piece NumberOfPoints=\""<<(int)wave_pos.size()<<"\" NumberOfCells=\""<<(int)wave_pos.size()-1<<"\">"<<endl;
  myfile<<endl;
  
  //COORDS. INFO.
  myfile<<"      <!-- Point coordinates -->"<<endl;
  myfile<<"      <Points>"<<endl;
  myfile<<"        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">"<<endl;
  for (int n=0;n<(int)wave_pos.size();n++)
    myfile<<"        "<<wave_pos[n]<<" "<<times[n]<<" 0.0"<<endl;    

  myfile<<"        </DataArray>"<<endl;
  myfile<<"      </Points>"<<endl;
  myfile<<endl;

  int count = 0; //for ensuring 6 columns of data per data type (e.g. coords)
  //CELLS (CONNECTIVITIES) INFO.
  myfile<<"      <!-- Connectivity (lines between points) -->"<<endl;
  myfile<<"      <Cells>"<<endl;
  myfile<<"        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">"<<endl;
  for (int n=0;n<(int)wave_pos.size()-1;n++)
    myfile<<"        "<<n<<" "<<n+1<<endl;

  myfile<<"        </DataArray>"<<endl;
  myfile<<"        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">"<<endl;
  myfile<<"        ";
  for (int n=1;n<(int)wave_pos.size();n++){
    count++;
    myfile<<2*n<<" ";    
    if (count % 6 == 0){
      myfile<<endl;
      myfile<<"        ";
    }
  }
  count = 0;
  myfile<<endl;
  myfile<<"        </DataArray>"<<endl;
  myfile<<"        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">"<<endl;
  myfile<<"        ";
  for (int n=0;n<(int)wave_pos.size();n++){
    count++;
    myfile<<3<<" ";
    if (count % 6 == 0){
      myfile<<endl;
      myfile<<"        ";
    }
  }
  count = 0;
  myfile<<endl;
  myfile<<"        </DataArray>"<<endl;
  myfile<<"      </Cells>"<<endl;

  //END OF FILE
  myfile<<"    </Piece>"<<endl;
  myfile<<"  </UnstructuredGrid>"<<endl;
  myfile<<"</VTKFile>"<<endl;
  

  return;
}
//-------------------------------------------------
void ParaviewWriter::WriteLinePlot(const char* &filename,vector<double> &xcoords,vector<double> &density,vector<double> &velocity,vector<double> &pressure){

  //ideal format: titleA titleB ... titleZ
  //              dataA  dataB  ... dataZ
  //              end of file
  
  ofstream myfile(filename);

  if (!myfile){ //checking if file opened successfully
    cerr<<"Error: Could Not Open File "<<filename<<endl;
    return;
  }

  //Column vector titles
  myfile<<"Analytical_Coords"<<"	"<<"Analytical_Density"<<"	"<<"Analytical_Velocity"<<"	"<<"Analytical_Pressure"<<endl;

  //Data
  for (int n=0;n<(int)xcoords.size();n++){
    myfile<<xcoords[n]<<"	"<<density[n]<<"	"<<velocity[n]<<"	"<<pressure[n]<<endl;

  }


  return;
}
//-------------------------------------------------
ParaviewWriter::~ParaviewWriter(){}
//-------------------------------------------------
