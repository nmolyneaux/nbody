/* ------------------------------------------------

Author: Nicholas Molyneaux
Date : 12 November 2015

Brute force approach to the n-body problem. Implemented
or the course of program parallelization for PC clusters

 ------------------------------------------------*/

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include <ctime>
#include <mpi.h>
#include <vector>
#include <stdexcept> // deal with exceptions

// --------------------------------------- 
// -------- Read into Eigen matrix -------
std::vector<double> readDataFile(const char *fileName)
{
  int nbCols = 0;
  int nbRows = 0;

  std::ifstream file;
  file.open(fileName);
  if (!file.is_open())
      throw std::invalid_argument("File not open !");

  // buffer for temp storage of values
  std::vector<double> buffer(0);
  std::string line;

  int tempCols = 0;

  while (std::getline(file, line))
    {      
      tempCols = 0;
      std::stringstream stream(line);
      double val;
      while( stream  >> val )
	{
	  buffer.push_back(val);
	  tempCols++;	  	  
	}
      if (nbCols == 0)
	nbCols = tempCols;
      nbRows++;
    }

  file.close();
  return buffer;
};

// ----------------------------------------
// ----------------- MAIN -----------------
int main(int argc, char* argv[])
{
  
  // --------------- Initialization ----------------
  // constant variables required by the model
  const double tol = 1e3;
  const double G = 6.674e-11;

  // Loads the data for the simulation into matrix. Separator is a comma ' ' (space)
  // Only the main processor does this, then sends the split data to workers
  int nbBodies;
  std::vector<double> data_read;
  std::vector<double> mass;
  std::vector<double> positions;
  std::vector<double> velocities;
  
  // --------------- Loading data  ----------------
  
  std::cout << "Loading data... ";    
  data_read = readDataFile(argv[1]);     
  nbBodies = data_read.size() / 5;
  if (nbBodies < 0)
      throw std::invalid_argument("Number of bodies smaller than number of nodes");
  for (int i = 0; i<nbBodies; i++)
  {
      mass.push_back(data_read[i*5]);
      positions.push_back(data_read[i*5+1]);
      positions.push_back(data_read[i*5+2]);
      velocities.push_back(data_read[i*5+3]);
      velocities.push_back(data_read[i*5+4]);
  }
  std::vector<double> positions_fixed(positions);
  
  std::cout << "done" << std::endl;
  
  // --------------- Dumper ----------------
  // Creates file to write data to
  std::string file_name = "output.csv"; 
  std::ofstream outputFile(file_name.c_str());
 
  outputFile.precision(10);
  outputFile << "t,px,py" << std::endl;
  for (int k = 0; k < nbBodies;k++) 
      outputFile << 0 << "," << positions[2*k] << "," <<  positions[2*k+1] << std::endl;
  
  // --------------- Time iterations ----------------
  double dt = 24*60*60;
  double time_max = 365*24*60*60;  
  double t = 0;
  
  
  for (; t < time_max; t += dt)
  {            
      //if (t % (10*24*60*60) == 0)
      std:: cout << "At time: " << t+dt << std::endl;
	   
      for(int i = 0; i < nbBodies; i++)
      {	  
	  for(int j = 0; j<nbBodies; j++)
	  {	      
	      double distance = sqrt(pow(positions_fixed[2*j]-positions[2*i],2) + pow(positions_fixed[2*j+1]-positions[2*i+1],2));
	      if (distance > tol)
	      {
		  velocities[2*i] += dt * (positions_fixed[2*j] - positions[2*i]) * G*mass[j]/(distance*distance*distance);
		  velocities[2*i+1] += dt * (positions_fixed[2*j+1] - positions[2*i+1]) * G*mass[j]/(distance*distance*distance);
	      }
	  }
	  
	  positions[2*i] += dt * velocities[2*i];
	  positions[2*i+1] += dt * velocities[2*i+1];
	  
	  
      }
      positions_fixed = positions;
      
      for (int k = 0; k < nbBodies;k++) 
	  outputFile << t+dt << "," << positions[2*k] << "," <<  positions[2*k+1] << std::endl;	 
  }
  
      // --------------- Finalization  ----------------
  outputFile.close();
  
}

