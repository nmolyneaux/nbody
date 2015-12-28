clear all; close all; clc;

N = 100000;
A=[];

    mass = [rand(N,1)*1e9 + 5e9];
    positions = rand(N,2) * 200 - 100; 
    velocities = rand(N,2) * 20;
    A=[A;mass,positions,velocities];
    
size(A)
dlmwrite('bodies_100000.dat',A,'delimiter',' ','precision',8);