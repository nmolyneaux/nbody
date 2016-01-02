clear all; close all; clc;delete(gcp('nocreate'));

N = 100000;
A = zeros(N,5);
pool = parpool(8);
parfor i=1:N
    mass = rand(1,1)*5e15 + 5e16;
    r = 325*(0.05 + rand(1,1));
    t = 2 * pi * rand(1,1);
    positions = [r * cos(t) ,r * sin(t)]; 
    velocities =  2*[-positions(:,2),positions(:,1)];    
    A(i,:) = [mass,positions,velocities];
end

size(A)
dlmwrite('bodies_100000.dat',A,'delimiter',' ','precision',8);