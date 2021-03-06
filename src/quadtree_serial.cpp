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
    std::clock_t total_time;
    std::clock_t start = std::clock();
    
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
    for (int i = 0; i < nbBodies; i++)
    {
	mass.push_back(data_read[i*5]);
	positions.push_back(data_read[i*5+1]);
	positions.push_back(data_read[i*5+2]);
	velocities.push_back(data_read[i*5+3]);
	velocities.push_back(data_read[i*5+4]);
    }
    std::vector<double> positions_fixed(positions);    
    std::cout << "done, " << nbBodies << " bodies." << std::endl;
    
    double dt = 0.001;
    double time_max = 10;
    double t = 0;   
    
    // --------------- Building first QT  ----------------
    Quadtree quad_tree =  Quadtree(0,0,10e11, 10e11, dt, 1, 0.5);    
    for (int i = 0; i < nbBodies; i++)
    {
	Body body(positions[2*i], positions[2*i+1], mass[i], velocities[2*i], velocities[2*i+1]);
	quad_tree.insertBody(body, quad_tree.root);
    }
    std::ofstream outputFile;
    /*
    outputFile.open(argv[2]);
    outputFile.precision(10);
    outputFile << "t,px,py" << std::endl;
    quad_tree.printPositions(quad_tree.root, t, outputFile);
    */
    // --------------- Time iterations  ----------------
    std::vector<double> bodies_data;
    std::cout << "Starting time loop" << std::endl;
    for (; t < time_max; t += dt)
      { 
	quad_tree.calculateForcesInBranch(quad_tree.root);	
	quad_tree.moveBodies(quad_tree.root);
	//quad_tree.printPositions(quad_tree.root, t+dt, outputFile);
	quad_tree.collectBodies(bodies_data, quad_tree.root, false);
        quad_tree.empty();

	// --------------- Re-building QT  ----------------
	for (int i = 0; i < nbBodies; i++)
	{
	    Body body(bodies_data[i*5], bodies_data[i*5+1], bodies_data[i*5+2], bodies_data[i*5+3], bodies_data[i*5+4]);
	    quad_tree.insertBody(body, quad_tree.root);
	}
	bodies_data.clear();
    }
    
    // --------------- Finalization  ----------------
    total_time = (std::clock() - start) / (double) CLOCKS_PER_SEC;
    std::cout << "time taken: " << total_time  << std::endl;
    outputFile.close();
}
