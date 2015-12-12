clear all; close all; clc;

N = 10000;

mass = zeros(N,1) + 1e9;%rand(N,1)*1e27 + 5e28;
positions = rand(N,2) * 1000;
velocities = rand(N,2) * 2;% 60000.;
A = [mass,positions,velocities];
%A = [A;1e32,0,0,0,0];
dlmwrite('bodies_10000.dat',A,'delimiter',' ');