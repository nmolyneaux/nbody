clear; close all; clc;

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
time = time(1:8000,:)
figure();
hold on;
plot(transpose(1:(size(time,1)-1)),time(1:end-1,1),'Color',[0    0.4470    0.7410])
plot(transpose(1:(size(time,1)-1)),time(1:end-1,2),'Color',[0.8500    0.3250    0.0980])
plot(transpose(1:(size(time,1)-1)),time(1:end-1,3),'Color',[0.9290    0.6940    0.1250])
plot(transpose(1:(size(time,1)-1)),time(1:end-1,4),'Color',[0.4940    0.1840    0.5560])
hold off;
