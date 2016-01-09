clear; close all; clc;

nodes = [1,2,4,8,12,16,20,24,28,32];

base_name = {'../results/bf-1000000_timing_nodes_'}
base_name_end = {'.csv_rank_'};

computeTimes = [];

for i=1:length(nodes)
    time = [];
    for j = 0:(nodes(i)-1)
        filename = strcat(base_name,num2str(nodes(i)),base_name_end,num2str(j),'.csv');
        data = dlmread(char(filename),',',0,1);
        time = [time,data];
        %if (j == 0 )
        %    disp('rank 0');
        %    data(end)
        %end
    end    
    computeTimes = [computeTimes,nodes(i).*time(3,:)];
    nodes(i)
    mean(time,2)
    std(time,[],2)
    time
end
computeTimes = computeTimes./44875;
%%

nodesCS = cumsum(nodes);
ymin = 0.75; ymax = 1.15;
figure()
hold on;
for i = 1:length(nodes)
    plot([nodesCS(i),nodesCS(i)],[ymin,ymax],'k')
end
scatter(1:147,computeTimes);
ax = gca;
ax.XTick = [0.5,2,5,11,21,35,53,75,101,131];
ax.XTickLabel = {'1','2','4','8','12','16','20','24','28','32'};
hold off;
dlmwrite('../figures/computeTimesNormalizedBF10e6.csv',[1:147;computeTimes]')




