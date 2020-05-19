%% SC42100 Assignment 1, exercise 5
% Maxime Croft (4390024) and Jessie van Dam (4395832)
clear all; close all; clc;

% Load given data
load -ascii traffic.mat
load -ascii flow.mat
load -ascii traveltime.mat
load -ascii capacities.mat

% Specifying graph edges (s,t) in node pairs
s = []; % tail of head
t = []; % tail of edge

% For 1, node is tail and thus in s
% For -1, node is head and thus in t
for i = 1:size(traffic,2)
    for j = 1:size(traffic, 1)
        if traffic(j,i) == 1
            s = [s; j];
        end
        if traffic(j,i) == -1
            t = [t; j];
        end
    end      
end

G = digraph(s,t);  % defining directed graph G without weights
W = adjacency(G,traveltime);

[dist,path,pred] = graphshortestpath(W,1,13);

%% Assignment A - shortest path between node 1 and 13

