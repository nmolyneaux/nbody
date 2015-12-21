clear all; close all; clc;

N = 100;
A=[];

    mass = [rand(100,1)*1e9 + 5e9];
    positions = [rand(100,1) * 20 + 50,rand(100,1) * 20 - 10]; 
    %velocities = [velocities;rand(1000,2) * 10-5];% 60000.;
    velocities =  2*[-positions(:,2),positions(:,1)-60];    
    A=[A;mass,positions,velocities];
    
    mass = [rand(100,1)*1e10 + 5e10];
    positions = 2+[rand(100,1) * 20 - 70,rand(100,1) * 20 - 10]; 
    %velocities = [velocities;rand(1000,2) * 10-5];% 60000.;
    velocities =  2*[-positions(:,2),positions(:,1)+60];    
    A=[A;mass,positions,velocities];
 
A = [A;1e15,60,0,0,0];
A = [A;1e15,-60,0,0,0];

size(A)
dlmwrite('bodies_202.dat',A,'delimiter',' ','precision',8);