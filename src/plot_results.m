clear all; close all; clc;

%% Loads data from csv written from nbdoy
data = dlmread('qt_serial.csv',',',1,0);

nbBodies = length(data(:,1))/length(unique(data(:,1)));

xlimits = max(abs(data(:,2)));
ylimits = max(abs(data(:,3)));
plotLimits = 1.1 * max(xlimits,ylimits);    


%% Plots the data for all time steps
for i = 1:length(unique(data(:,1)))
    plot(data(1+nbBodies*(i-1):nbBodies+nbBodies*(i-1),2),data(1+nbBodies*(i-1):nbBodies+nbBodies*(i-1),3),'o');
    xlim([-plotLimits,plotLimits]);ylim([-plotLimits,plotLimits]);daspect([1,1,1]);
    xlim([-10e11,10e11]);ylim([-10e11,10e11]);daspect([1,1,1]);    
    pause(0.01);
end

