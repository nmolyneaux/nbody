clear all; close all; clc;

N = 1500;
N1 = N;
A=[];

for i=1:N
    mass = rand(1,1)*1e11 + 1e11;
    r = 12 + 1*(0 + rand(1,1));
    t = 2 * pi * rand(1,1);
    positions = [r * cos(t) ,r * sin(t)]; 
    velocities =  15*[-positions(:,2),positions(:,1)];  
    velocities(:,1) = velocities(:,1);% + 250;
    A=[A;mass,positions,velocities];
end

for i=1:N
    mass = rand(1,1)*1e11 + 1e11;
    r = 12 + 1*(0 + rand(1,1));
    t = 2 * pi * rand(1,1);
    positions = [r * cos(t) ,r * sin(t)]; 
    velocities =  15*[-positions(:,2) ,positions(:,1)];    
    velocities(:,1) = velocities(:,1) - 250;
    positions(:,1) = positions(:,1) + 50; 
    A=[A;mass,positions,velocities];
end

A = [A;2e17,0,0,0,0];
A = [A;2e17,50,0,-250,0];

size(A)
dlmwrite('bodies_10000.dat',A,'delimiter',' ','precision',8);