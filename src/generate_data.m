clear all; close all; clc;

N = 1001;
N1 = N;
A=[];

for i=1:N
    mass = rand(1,1)*5e15 + 5e16;
    r = 325*(0.05 + rand(1,1));
    t = 2 * pi * rand(1,1);
    positions = [r * cos(t) ,r * sin(t)]; 
    velocities =  2*[-positions(:,2),positions(:,1)];    
    A=[A;mass,positions,velocities];
end


A = [A;6e17,0,0,0,0];
%A = [A;1e20,-250,0,0,0];

size(A)
dlmwrite('bodies_1000.dat',A,'delimiter',' ','precision',8);