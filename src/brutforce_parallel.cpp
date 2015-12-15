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
  MPI_Init(&argc, &argv); 

  double start_time, end_time,start_time_total;
  start_time_total = MPI_Wtime();
  start_time = MPI_Wtime();

  
  int my_rank;
  int nb_proc;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); 
  MPI_Comm_size(MPI_COMM_WORLD, &nb_proc); 
  
  std::ofstream timingFile(argv[3]);
  if (my_rank == 0)
  {
      std::ofstream timingFile(argv[3]);
      timingFile.precision(10);
  }

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
  if (my_rank == 0)
  {           
      std::cout << "Loading data... ";    
      data_read = readDataFile(argv[1]);     
      nbBodies = data_read.size() / 5;
      if (nbBodies < nb_proc)
	  throw std::invalid_argument("Number of bodies smaller than number of nodes");
      for (int i = 0; i<nbBodies; i++)
      {
	  mass.push_back(data_read[i*5]);
	  positions.push_back(data_read[i*5+1]);
	  positions.push_back(data_read[i*5+2]);
	  velocities.push_back(data_read[i*5+3]);
	  velocities.push_back(data_read[i*5+4]);	  
      }
      std::cout << "done" << std::endl;
  }
  // Broadcasts the number of bodies to all nodes
  MPI_Bcast(&nbBodies, 1, MPI_INT, 0, MPI_COMM_WORLD);
  
  // --------------- Sending initial data ----------------
  // calculating block length for all workers
  int *block_length = new int[nb_proc];
  int *start_position = new int[nb_proc];
  for (int i = 0; i < (nb_proc-1); i++)
  {
      block_length[i] = 2 * (nbBodies/nb_proc);
      start_position[i] = i * 2 * (nbBodies/nb_proc);
  }
  start_position[nb_proc-1] = (nb_proc-1) * 2 * (nbBodies/nb_proc);
  block_length[nb_proc-1] = 2 * (nbBodies - (nbBodies/nb_proc)*(nb_proc-1));

  // Sending initial data to all workers (including 0).
  // For all workers resize the vectors to contain the data 
  if (my_rank != 0)
  {
      mass.resize(nbBodies);
      positions.resize(2*nbBodies);
  }
  
  // Local vectors to store the positions and velocities of the "local" bodies on each node.
  // once the positions_local vector is updated, it will be sent to all nodes
  std::vector<double> velocities_local(block_length[my_rank]);   
  std::vector<double> positions_local(block_length[my_rank]);   
  
  end_time = MPI_Wtime();
  if (my_rank == 0)
      timingFile << "loading_data, " << MPI_Wtime() - start_time << std::endl;

  // Actually sends data to all nodes
  start_time = MPI_Wtime();
  MPI_Bcast(&mass[0], nbBodies, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&positions[0], 2*nbBodies, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatter(&positions[0], block_length[my_rank], MPI_DOUBLE, &positions_local[0], block_length[my_rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatter(&velocities[0], block_length[my_rank], MPI_DOUBLE, &velocities_local[0], block_length[my_rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

  end_time = MPI_Wtime();
  if (my_rank == 0)   
      timingFile << "sending_initial," << end_time - start_time << std::endl;

  // --------------- Dumper ----------------
  // Creates file to write data to
  start_time = MPI_Wtime();
  std::string file_name = argv[2]; 
  std::ofstream outputFile(file_name.c_str());
  if (my_rank == 0)
  {
      std::ofstream outputFile(file_name.c_str());
      outputFile.precision(10);
      outputFile << "t,px,py" << std::endl;
      for (int k = 0; k < nbBodies;k++) 
	  outputFile << 0 << "," << positions[2*k] << "," <<  positions[2*k+1] << std::endl;
  }
  
  // --------------- Time iterations ----------------
  double dt = 0.1;
  double time_max = 0.1;  
  double t = 0;
  
  end_time = MPI_Wtime();
  if (my_rank == 0)   
  {
      timingFile << "writing_inital," << end_time - start_time << std::endl;  
      timingFile << "forces,allgatherv,writing" << std::endl;
  }

  for (; t < time_max; t += dt)
  {            
      if (my_rank == 0)
      {
	  std:: cout << "At time: " << t << std::endl;	   
      }
      start_time = MPI_Wtime();

      for(int i = 0; i < 0.5 * block_length[my_rank]; i++)
      {	  
	  for(int j = 0; j<nbBodies; j++)
	  {	      
	      double distance = sqrt(pow(positions[2*j]-positions_local[2*i],2) + pow(positions[2*j+1]-positions_local[2*i+1],2));
	      if (distance > tol)
	      {
		  velocities_local[2*i] += dt * (positions[2*j] - positions_local[2*i]) * G*mass[j]/(distance*distance*distance);
		  velocities_local[2*i+1] += dt * (positions[2*j+1] - positions_local[2*i+1]) * G*mass[j]/(distance*distance*distance);
	      }
	  }
	  positions_local[2*i] += dt * velocities_local[2*i];
	  positions_local[2*i+1] += dt * velocities_local[2*i+1];
      }
      
      end_time = MPI_Wtime();
      if (my_rank == 0)   
	  timingFile << end_time - start_time << ",";

      start_time = MPI_Wtime();

      MPI_Allgatherv(&positions_local[0], block_length[my_rank], MPI_DOUBLE, &positions[0], block_length, start_position,  MPI_DOUBLE, MPI_COMM_WORLD);

      end_time = MPI_Wtime();
      if (my_rank == 0)   
	  timingFile << end_time - start_time << ",";
      
      start_time = MPI_Wtime();
      if (my_rank == 0)
      {	  
	  for (int k = 0; k < nbBodies;k++) 
	      outputFile << t+dt << "," << positions[2*k] << "," <<  positions[2*k+1] << std::endl;	 
      }
      end_time = MPI_Wtime();
      if (my_rank == 0)   
	  timingFile << end_time - start_time << std::endl;
      
  }
  
  // --------------- Finalization  ----------------
  if (my_rank == 0)
  {
      outputFile.close();
      timingFile << "total_time," << MPI_Wtime() - start_time_total << std::endl;
      timingFile.close();
  }
  
  delete[] start_position;
  delete[] block_length;

  MPI_Finalize();
}
