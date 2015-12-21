clear all; close all; clc;

N = 1000;
A=[];

    mass = [rand(N,1)*1e9 + 5e9];
    positions = rand(N,2) * 200 - 100; 
    velocities = rand(N,2) * 20;
    A=[A;mass,positions,velocities];
    
size(A)
dlmwrite('bodies_1000.dat',A,'delimiter',' ','precision',8);