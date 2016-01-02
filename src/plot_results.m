clear all; clc; close all;

base_name = {'../results/qt-two-masses-timing_rank_'}


time = [];
for j = 0:3
    filename = strcat(base_name,num2str(j),'.csv')
    data = dlmread(char(filename),',',0,1);
    time = [time,data];
    if (j == 0 )
        disp('rank 0');
        data(end)
    end
end
time = time(1:8000,:);
figure();
hold on;
plot(transpose(1:(size(time,1)-1)),time(1:end-1,1),'Color',[0    0.4470    0.7410])
plot(transpose(1:(size(time,1)-1)),time(1:end-1,2),'Color',[0.8500    0.3250    0.0980])
plot(transpose(1:(size(time,1)-1)),time(1:end-1,3),'Color',[0.9290    0.6940    0.1250])
plot(transpose(1:(size(time,1)-1)),time(1:end-1,4),'Color',[0.4940    0.1840    0.5560])
hold off;


%% Loads data from csv written from nbdoy
data = dlmread('../results/qt-two-masses-results.csv',',',1,0);
%%
%nbBodies = length(data(:,1))/length(unique(data(:,1)));
nbBodies = 3002;
%%
xlimits = max(abs(data(:,2)));
ylimits = max(abs(data(:,3)));
plotLimits =  0.1*max(xlimits,ylimits);

nbProcs = length(unique(data(:,5)));
% Plots the data for all time steps

%vidObj = VideoWriter('peaks.avi');
%vidObj.FrameRate = 10;
%vidObj.VideoCompressionMethod = 'MPEG-4';
%vidObj.Height = 1080;
%vidObj.Width = 1920;

%open(vidObj);   

%h = figure('units','pixels','position',[0,0 1920 1080]);
figure();
%%
bodies_per_node=[];
colors = zeros(nbBodies,3);
for i = 2:length(unique(data(:,1)))-1
    %hold on
    %subplot(2,1,1,'position', [0.05, 0.3, 0.9, 0.7] )
    %set(gca,'nextplot','replacechildren');
    [a,~] = hist(data(1+nbBodies*(i-1):nbBodies+nbBodies*(i-1),5),unique(data(1+nbBodies*(i-1):nbBodies+nbBodies*(i-1),5)));
    bodies_per_node = [bodies_per_node;a];
    
    colors(data(1+nbBodies*(i-1):nbBodies+nbBodies*(i-1),5)==0,:) = repmat([0    0.4470    0.7410],a(1),1);
    colors(data(1+nbBodies*(i-1):nbBodies+nbBodies*(i-1),5)==1,:) = repmat([0.8500    0.3250    0.0980],a(2),1);
    colors(data(1+nbBodies*(i-1):nbBodies+nbBodies*(i-1),5)==2,:) = repmat([0.9290    0.6940    0.1250],a(3),1);
    colors(data(1+nbBodies*(i-1):nbBodies+nbBodies*(i-1),5)==3,:) = repmat([0.4940    0.1840    0.5560],a(4),1);

    scatter(data(1+nbBodies*(i-1):nbBodies+nbBodies*(i-1),2),data(1+nbBodies*(i-1):nbBodies+nbBodies*(i-1),3),[],colors,'filled');
    %scatter(data(1+nbBodies*(i-1):nbBodies+nbBodies*(i-1),2),data(1+nbBodies*(i-1):nbBodies+nbBodies*(i-1),3));
    xlim([-200,700]);ylim([-200,200]);daspect([1,1,1]);
    title(strcat('Time step number ',num2str(i)))
    %subplot(2,1,2,'position', [0.05, 0.05, 0.9, 0.2] )
    %histogram(data(1+nbBodies*(i-1):nbBodies+nbBodies*(i-1),5),'BinMethod', 'integers')
    %ylim([0,1.2*nbBodies/nbProcs]);
    pause(0.00001)
    %hold off
    %currFrame = getframe(gcf);
    %writeVideo(vidObj,currFrame);
end
%%
figure();
plot(bodies_per_node);
ylim([0,800]);

%close(vidObj);


