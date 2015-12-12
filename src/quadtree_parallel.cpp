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
#include <numeric> // cumsum
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

  double start_timer = MPI_Wtime();

  // Loads the data for the simulation into matrix. Separator is a comma ' ' (space)
  // Only the main processor does this, then sends the split data to workers
  int nbBodies;
  std::vector<double> bodies_data;
  std::vector<double> mass;
  std::vector<double> positions;
  std::vector<double> velocities;

  // --------------- Loading data  ----------------
  if (my_rank == 0)
  {           
      std::cout << "Loading data... ";    
      bodies_data = readDataFile(argv[1]);     
      nbBodies = bodies_data.size() / 5;
      if (nbBodies < nb_proc)
	  throw std::invalid_argument("Number of bodies smaller than number of nodes");
      for (int i = 0; i<nbBodies; i++)
      {
	  mass.push_back(bodies_data[i*5]);
	  positions.push_back(bodies_data[i*5+1]);
	  positions.push_back(bodies_data[i*5+2]);
	  velocities.push_back(bodies_data[i*5+3]);
	  velocities.push_back(bodies_data[i*5+4]);	  
      }
      std::cout << "done" << std::endl;
  }
  // Broadcasts the number of bodies to all nodes
  MPI_Bcast(&nbBodies, 1, MPI_INT, 0, MPI_COMM_WORLD);
  
  bodies_data.resize(nbBodies*5);
  MPI_Bcast(&bodies_data[0], 5*nbBodies, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  
  double dt = 0.1;
  double time_max = 0.1;
  double t = 0;   
  
  Quadtree quad_tree =  Quadtree(0,0,10e4, 10e4, dt, nb_proc);
    
  for (int i = 0; i < nbBodies; i++)
  {
      Body body(bodies_data[i*5+1], bodies_data[i*5+1+2], bodies_data[i*5], bodies_data[i*5+3], bodies_data[i*5+4]);     
      quad_tree.insertBody(body, quad_tree.root);	  
  }
  
  std::vector<std::vector < Node *> > node_assignment(nb_proc);
  
  std::ofstream outputFile("qt_parallel.csv");
  if (my_rank==0)
  {
      outputFile.precision(10);
      outputFile << "t,px,py" << std::endl;
      quad_tree.printPositions(quad_tree.root, t, outputFile);//
      std::cout << "Starting time loop" << std::endl;
  }
  
  int *nb_bodies_per_node = new int[nb_proc];
  int *block_length = new int[nb_proc];
  int *start_position = new int[nb_proc];

  std::vector<double> bodies_data_local;
  for (; t < time_max; t += dt)
  {   	
      if (my_rank==0)
	  std::cout << "We are at time: " << t << std::endl;

      node_assignment = quad_tree.findLocalNodes(nb_bodies_per_node);

      std::cout << node_assignment[my_rank].size() << std::endl;

      start_position[0] = 0;
      block_length[0] = (nb_bodies_per_node[0])*5;
      for (int i = 0; i < nb_proc -1; i++)
      {
	  block_length[i] = (nb_bodies_per_node[i])*5;
	  start_position[i+1] = (nb_bodies_per_node[i])*5;
      }

      for (int i = 0; i < node_assignment[my_rank].size(); i++)	  
      {
	  quad_tree.calculateForcesInBranch( (*node_assignment[my_rank][i]) );
      }      
      
      bodies_data_local.clear();
      for (int i = 0; i < node_assignment[my_rank].size(); i++)	  
      {
	  quad_tree.moveBodies( (*node_assignment[my_rank][i]) );
	  quad_tree.collectBodies(bodies_data_local, (*node_assignment[my_rank][i]) );
      }     

      //bodies_data.clear();
      std::cout << block_length[my_rank] << " " << block_length[0] << " " << start_position[0] << std::endl;
      MPI_Allgatherv(&bodies_data_local[0], block_length[my_rank], MPI_DOUBLE, &bodies_data[0], block_length, start_position,  MPI_DOUBLE, MPI_COMM_WORLD);
      
      if(my_rank == 0)
	  quad_tree.printPositions(quad_tree.root, t+dt, outputFile);
      
      quad_tree.empty();
      
      std::cout << bodies_data.size() << std::endl;
      for (int j = 0; j < bodies_data.size(); j++)
	  std::cout << bodies_data[j] << " ";
      std::cout << std::endl;

      for (int i = 0; i < nbBodies; i++)
      {
	  Body body(bodies_data[i*5], bodies_data[i*5+1], bodies_data[i*5+2], bodies_data[i*5+3], bodies_data[i*5+4]);
	  std::cout << "inserting " << body << std::endl;
	  quad_tree.insertBody(body, quad_tree.root);
	  }
  }
  MPI_Finalize();
  
}


