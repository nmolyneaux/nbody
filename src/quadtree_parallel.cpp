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
  int my_rank;
  int nb_proc;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); 
  MPI_Comm_size(MPI_COMM_WORLD, &nb_proc); 

  std::ofstream timingFile(argv[3]);
  timingFile.precision(10);

  double start_time, end_time;
  start_time = MPI_Wtime();
  // Loads the data for the simulation into matrix. Separator is a comma ' ' (space)
  // Only the main processor does this, then sends the split data to workers
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
  timingFile << "loading_data, " << MPI_Wtime() - start_time << std::endl;

  // Broadcasts the number of bodies to all nodes
  MPI_Bcast(&nbBodies, 1, MPI_INT, 0, MPI_COMM_WORLD);
  bodies_data.resize(nbBodies*5);
  
  start_time = MPI_Wtime();
  MPI_Bcast(&bodies_data[0], 5*nbBodies, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  end_time = MPI_Wtime();
  if (my_rank == 0)      
      timingFile << "sending_initial," << end_time - start_time << std::endl;
  
  double dt = 0.1;
  double time_max = 0.5;
  double t = 0;   
  
  start_time = MPI_Wtime();
  Quadtree quad_tree =  Quadtree(0,0,10e12, 10e12, dt, nb_proc);
  for (int i = 0; i < nbBodies; i++)
  {
      Body body(bodies_data[i*5+1], bodies_data[i*5+2], bodies_data[i*5], bodies_data[i*5+3], bodies_data[i*5+4]);     
            quad_tree.insertBody(body, quad_tree.root);
  }
  end_time = MPI_Wtime();
  if (my_rank == 0)
      timingFile << "builing_first_qt," << end_time - start_time << std::endl;

  std::vector<std::vector < Node *> > node_assignment(nb_proc);
  
  std::ofstream outputFile(argv[2]);
  if (my_rank==0)
  {
      timingFile << "forces,allgatherv,writing,rebuilding" << std::endl;
      outputFile.precision(10);
      outputFile << "t,px,py" << std::endl;
      quad_tree.printPositions(quad_tree.root, t, outputFile);
      std::cout << "Starting time loop" << std::endl;
  }
  

  int *nb_bodies_per_node = new int[nb_proc];
  int *block_length = new int[nb_proc];
  int *start_position = new int[nb_proc];

  std::vector<double> bodies_data_local;
  for (; t < time_max; t += dt)
  {   	
      if (my_rank==0)
      
	  std::cout << "We are at time: " << t+dt << std::endl;

      // Load balancing and arrays containg the sizes for MPI gatehring
      start_time = MPI_Wtime();

      node_assignment = quad_tree.findLocalNodes(nb_bodies_per_node);
      start_position[0] = 0;
      //block_length[0] = (nb_bodies_per_node[0])*5;
      for (int i = 0; i < (nb_proc - 1); i++)
      {
	  *(block_length + i) = *(nb_bodies_per_node + i)*5;
	  *(start_position + i+1) = *(start_position + i) + *(nb_bodies_per_node + i)*5;
      }
      *(block_length + (nb_proc -1)) = *(nb_bodies_per_node + (nb_proc -1))*5;
      
      // Calculates acceleration acting oin all bodies
      for (int i = 0; i < node_assignment[my_rank].size(); i++)	  
      {
	  quad_tree.calculateForcesInBranch( (*node_assignment[my_rank][i]) );
      }      
      
      // Move and then collect the new positions
      bodies_data_local.clear();
      for (int i = 0; i < node_assignment[my_rank].size(); i++)	  
      {
	  quad_tree.moveBodies( (*node_assignment[my_rank][i]) );
	  quad_tree.collectBodies(bodies_data_local, (*node_assignment[my_rank][i]) );
      }     

      // Give all nodes the new information about the positoins and velocities
      bodies_data.clear();

      end_time = MPI_Wtime();
      if (my_rank == 0)
	  timingFile << end_time - start_time;

      start_time = MPI_Wtime();
      MPI_Allgatherv(&bodies_data_local[0], block_length[my_rank], MPI_DOUBLE, &bodies_data[0], block_length, start_position,  MPI_DOUBLE, MPI_COMM_WORLD);
      end_time = MPI_Wtime();
      if (my_rank == 0)
	  timingFile << "," << end_time - start_time;

      start_time = MPI_Wtime();

      // print time setp to file
      if(my_rank == 0)
	  quad_tree.printPositions(quad_tree.root, t+dt, outputFile);
      end_time = MPI_Wtime();
      if (my_rank == 0)
	  timingFile << "," << end_time - start_time ;

      // rebuild tree for next time step
      start_time = MPI_Wtime();
      quad_tree.empty();        
      for (int i = 0; i < nbBodies; i++)
      {
	  Body body(bodies_data[i*5], bodies_data[i*5+1], bodies_data[i*5+2], bodies_data[i*5+3], bodies_data[i*5+4]);     
	  quad_tree.insertBody(body, quad_tree.root);
      }
      end_time = MPI_Wtime();
      if (my_rank == 0)
	  timingFile << "," << end_time - start_time << std::endl;

  }
  MPI_Finalize();  
}


