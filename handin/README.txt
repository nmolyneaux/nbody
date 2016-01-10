To compile the code on a linux system with a distribution of CMake, MPI and C++ installed, it should be sufficent to:

cd nbody/src/
cmake CMakeLists.txt
make  

CMake will search for the MPI package and use it to compile the executables.
This should compile all four executables, i.e.

bruteforce-serial
bruteforce-parallel
quadtree-serial
quadtree-parallel

To run the code, command line arguments must be passed in the following order:
1) the name (and path) to a data file. They are contained in the folder data
2) the name to give the output file (even if this is commented on compilatiom , it is suggested to give a dummy name)
3) ONLY *-parallel executables, the name to give the file to store the time taken for each operation.
