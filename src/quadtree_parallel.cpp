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
// -------- Read data from csv -------
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
        
    std::ofstream timingFile;
    if (my_rank == 0)
    {
	std::string timing_file_name = argv[3];    
        timingFile.open(timing_file_name.c_str());
	timingFile.precision(10);
    }

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
    
    double start_time, end_time, start_time_total;
    start_time_total = MPI_Wtime();

    // Broadcasts the number of bodies to all nodes
    MPI_Bcast(&nbBodies, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    bodies_data.resize(nbBodies * 5);
    MPI_Bcast(&bodies_data[0], 5*nbBodies, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (my_rank == 0) timingFile << "sending_initial," << MPI_Wtime() - start_time << std::endl;
    
    double dt = 0.00005;
    double time_max = 0.3;
    double t = 0;   
    
    // --------------- Building first QT  ----------------
    Quadtree quad_tree =  Quadtree(0,0,10e11, 10e11, dt, nb_proc, 0.5);
    for (int i = 0; i < nbBodies; i++)
    {
	Body body(bodies_data[i*5+1], bodies_data[i*5+2], bodies_data[i*5], bodies_data[i*5+3], bodies_data[i*5+4]);     
	quad_tree.insertBody(body, quad_tree.root);
    }
    if (my_rank == 0) timingFile << "builing_first_qt," << MPI_Wtime() - start_time << std::endl;
    
    std::ofstream outputFile;
    
    if (my_rank==0)
    {      
        outputFile.open(argv[2]);
	outputFile.precision(10);
	outputFile << "t,px,py" << std::endl;
	quad_tree.printPositions(quad_tree.root, t, outputFile);	
    }
         
    
    int *nb_bodies_per_node = new int[nb_proc];
    int *block_length = new int[nb_proc];
    int *start_position = new int[nb_proc];
    
    bool with_proc = true;
    int doubles_per_body;

    if (with_proc)
	doubles_per_body = 6;
    else
	doubles_per_body = 5;
    
    bodies_data.resize(nbBodies * doubles_per_body);        
    std::vector<double> bodies_data_local;    
    int n = 0;
    
    // --------------- Starting time iterations  ----------------
    if (my_rank == 0) std::cout << "Starting time loop" << std::endl;
    
    for (; t < time_max; t += dt)
    {   	
	n += 1;
	if (my_rank==0 && (n%20) == 0)      
	    std::cout << "We are at time: " << t+dt << std::endl;
	
	// Load balancing and arrays containg the sizes for MPI gathering
	std::vector<std::vector < Node *> > node_assignment(nb_proc);
    
	node_assignment = quad_tree.findLocalNodes(nb_bodies_per_node);
	start_position[0] = 0;

	for (int i = 0; i < (nb_proc - 1); i++)
	{
	    *(block_length + i) = *(nb_bodies_per_node + i)*doubles_per_body;
	    *(start_position + i+1) = *(start_position + i) + *(nb_bodies_per_node + i)*doubles_per_body;
	}
	*(block_length + (nb_proc -1)) = *(nb_bodies_per_node + (nb_proc -1))*doubles_per_body;
	
	// Calculates acceleration acting on all bodies
	for (int i = 0; i < node_assignment[my_rank].size(); i++)	  
	{
	    quad_tree.calculateForcesInBranch( (*node_assignment[my_rank][i]), my_rank);
	}      
	
	// Move and then collect the new positions
	bodies_data_local.clear();
	for (int i = 0; i < node_assignment[my_rank].size(); i++)	  
	{
	    quad_tree.moveBodies( (*node_assignment[my_rank][i]) );
	    quad_tree.collectBodies(bodies_data_local, (*node_assignment[my_rank][i]), with_proc);
	}     
	
	// Give all nodes the new information about the positions and velocities
	bodies_data.clear();	
		
	MPI_Allgatherv(&bodies_data_local[0], block_length[my_rank], MPI_DOUBLE, &bodies_data[0], block_length, start_position,  MPI_DOUBLE, MPI_COMM_WORLD);
	
	// rebuild tree for next time step
	quad_tree.empty();
	//std::cout << quad_tree << std::endl;
	for (int i = 0; i < nbBodies; i++)
	{
	    //Body *body = new Body(bodies_data[i*doubles_per_body], bodies_data[i*doubles_per_body+1], bodies_data[i*doubles_per_body+2], bodies_data[i*doubles_per_body+3], bodies_data[i*doubles_per_body+4]);

	    Body body(bodies_data[i*doubles_per_body], bodies_data[i*doubles_per_body+1], bodies_data[i*doubles_per_body+2], bodies_data[i*doubles_per_body+3], bodies_data[i*doubles_per_body+4]);
	    quad_tree.insertBody(body, quad_tree.root);
	}
	/*
	// print time step to file
	
	if(my_rank == 0 && with_proc && (n%2) == 0)
	{
	    for (int i = 0; i < nbBodies; i++)
	    {
		outputFile << t+dt << "," << bodies_data[i*doubles_per_body] << "," << bodies_data[i*doubles_per_body+1] << "," << bodies_data[i*doubles_per_body+2] << "," << bodies_data[i*doubles_per_body+5] << std::endl;
	    }
	}
	else if (my_rank == 0 && !with_proc && (n%2) == 0)
	{
	    for (int i = 0; i < nbBodies; i++)
	    {
		outputFile <<  t+dt << "," << bodies_data[i*doubles_per_body] << "," << bodies_data[i*doubles_per_body+1] << "," << bodies_data[i*doubles_per_body+2] << std::endl;
	    }
	}
	*/
    }

    // --------------- Finalization  ----------------

    delete[] nb_bodies_per_node;
    delete[] block_length;
    delete[] start_position;
    
    if (my_rank==0)
    {      
	outputFile.close();
	timingFile << "total_time," << MPI_Wtime() - start_time_total << std::endl;
	timingFile.close();
    }
    MPI_Finalize();  
}
