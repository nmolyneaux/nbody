clear; close all; clc;

nodes = [1,2,4,8,12,16,20,24,28,32];

base_name = {'../results/bf-100000_timing_nodes_'}
base_name_end = {'.csv_rank_'};

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
    nodes(i)
    mean(time,2)
    std(time,[],2)
    min(time(4,:))
end