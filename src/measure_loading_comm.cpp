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
#include <iterator>
#include <stdexcept> // deal with exceptions

#include "Quadtree.hpp"
// --------------------------------------- 
// ------------ Read into buffer ---------
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
  MPI_Init(&argc, &argv); 

  double start_time, end_time,start_time_total;
  start_time_total = MPI_Wtime();
  start_time = MPI_Wtime();

  int my_rank;
  int nb_proc;

  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); 
  MPI_Comm_size(MPI_COMM_WORLD, &nb_proc); 
  
  std::ostringstream stream;
  stream << argv[2] << "_rank_" << my_rank << ".csv";
  std::string timing_file_name = stream.str();
  std::ofstream timingFile(timing_file_name.c_str());
   
  int nbBodies;
  std::vector<double> bodies_data;

  // --------------- Loading data  ----------------
  if (my_rank == 0)
  {           
      std::cout << "Loading data... ";    
      bodies_data = readDataFile(argv[1]);     
      nbBodies = bodies_data.size() / 5;
      if (nbBodies < nb_proc)
	  throw std::invalid_argument("Number of bodies smaller than number of nodes");
      std::cout << "done" << std::endl;
  }

  timingFile << "loading_data," << MPI_Wtime() - start_time << std::endl;
      
  // Broadcasts the number of bodies to all nodes
  start_time = MPI_Wtime();

  MPI_Bcast(&nbBodies, 1, MPI_INT, 0, MPI_COMM_WORLD); 
  bodies_data.resize(nbBodies * 5);
  MPI_Bcast(&bodies_data[0], 5*nbBodies, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  end_time = MPI_Wtime();
  timingFile << "sending_initial," << end_time - start_time << std::endl;
  
  
  timingFile << "total_time," << MPI_Wtime() - start_time_total << std::endl;
  timingFile.close();
  
  MPI_Finalize();  
}



