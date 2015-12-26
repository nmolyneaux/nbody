clear all; close all; clc;

%% Loads data from csv written from nbdoy
data = dlmread('qt_serial.csv',',',1,0);
%%
nbBodies = length(data(:,1))/length(unique(data(:,1)));
nbBodies = 1000;
%%
xlimits = max(abs(data(:,2)));
ylimits = max(abs(data(:,3)));
plotLimits =  0.1*max(xlimits,ylimits);

%nbProcs = length(unique(data(:,5)));
% Plots the data for all time steps

%vidObj = VideoWriter('peaks.avi');
%vidObj.FrameRate = 10;
%vidObj.VideoCompressionMethod = 'MPEG-4';
%vidObj.Height = 1080;
%vidObj.Width = 1920;

%open(vidObj);   

%h = figure('units','pixels','position',[0,0 1920 1080]);
figure();
for i = 1:length(unique(data(:,1)))
    %subplot(2,1,1)
    %set(gca,'nextplot','replacechildren');
    %subplot('position', [0.05, 0.3, 0.9, 0.7] );
    %scatter(data(1+nbBodies*(i-1):nbBodies+nbBodies*(i-1),2),data(1+nbBodies*(i-1):nbBodies+nbBodies*(i-1),3),[],data(1+nbBodies*(i-1):nbBodies+nbBodies*(i-1),5),'filled');
    scatter(data(1+nbBodies*(i-1):nbBodies+nbBodies*(i-1),2),data(1+nbBodies*(i-1):nbBodies+nbBodies*(i-1),3));

    xlim([-1050,1050]);ylim([-1050,1050]);daspect([1,1,1]);
    
    %subplot(2,1,2)
    %subplot('position', [0.05, 0.05, 0.9, 0.2] );
    %histogram(data(1+nbBodies*(i-1):nbBodies+nbBodies*(i-1),5),'BinMethod', 'integers')
    %ylim([0,1.2*nbBodies/nbProcs]);
    pause(0.01)
    %currFrame = getframe(gcf);
    %writeVideo(vidObj,currFrame);
end

%close(vidObj);


