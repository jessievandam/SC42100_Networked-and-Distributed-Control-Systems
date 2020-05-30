%% SC42100 Assignment 1, exercise 5
% Maxime Croft (4390024) and Jessie van Dam (4395832)
clear all; close all; clc;

% Load given data for traffic network
load -ascii traffic.mat     % node-link incidence matrix 
load -ascii flow.mat        % flow
load -ascii traveltime.mat  % capacities in vehicles per hour
load -ascii capacities.mat  % the minimum traveling times

% Specifying graph edges (s,t) in node pairs
s = []; % tail of head
t = []; % tail of edge

% For 1, node is tail and thus in s
% For -1, node is head and thus in t
for i = 1:size(traffic,2)
    for j = 1:size(traffic,1)
        if traffic(j,i) == 1
            s = [s; j];
        end
        if traffic(j,i) == -1
            t = [t; j];
        end
    end      
end

%% Assignment A - shortest path between node 1 and 13
% Adjacency matrix W with traveltime as weights
W = sparse(s, t, traveltime);
W = [W; zeros(1, size(W,2))]; % making W square

% Calculate shortest path from node 1 to 13
[dist,path,pred]  = graphshortestpath(W,1,13);

%% Assignment B - maximum flow between node 1 and 13
% Adjacency matrix W with capacities as weights
W = sparse(s, t, capacities);
W = [W; zeros(1, size(W,2))]; % making W square

% Calculate maximum flow from node 1 to 13
[MaxFlow, FlowMatrix, Cut] = graphmaxflow(W, 1, 13);

%% Assignment C - external inflow/outflow at each node
% Compute inflow (+) or outflow (-) of each node
B = traffic;
f = flow;
InOutFlow = B*f;

%% Assignment D - social optimum with respect to delays
% Defining variables
l = traveltime;
C = capacities;
M = size(B, 2);

% Minimization problem to compute social optimum
cvx_begin 
    variable f(M) 
    minimize sum((l.*C).*inv_pos(ones(M,1)-f./C)-l.*C)
    subject to 
        B(1,:)*f == 3329
        B(13,:)*f == -3329
        B(2:12,:)*f == 0
        0 <= f <= C
cvx_end

% Social optimum flow
f_SocialOptimum = f;

%% Assignment E - Wardrop equilibrium
% Minimization problem to compute Wardrop equilibrium
clear f;
cvx_begin
    variable f(M)   
    minimize sum(-l.*C.*log(1-f./C))
    subject to 
        B(1,:)*f == 3329
        B(13,:)*f == -3329
        B(2:12,:)*f == 0
        0 <= f <= C 
cvx_end 

% Wardrop Equilibrium flow
f_WardropEquilibrium = f;

%% Assignment F - Wardrop equilibrium with tolls
% Compute toll omega
omega = f_SocialOptimum.*(l./C)./((ones(M,1)-f_SocialOptimum./C).^2);

% Minimization problem to compute Wardrop equilibrium with tolls
clear f;

cvx_begin
    variable f(M)   
    minimize sum(-l.*C.*log(1-f./C)+f.*omega)
    subject to 
        B(1,:)*f == 3329
        B(13,:)*f == -3329
        B(2:12,:)*f == 0
        0 <= f <= C 
cvx_end 

% Wardrop Equilibrium with tollsflow
f_Tolls = f;

%% Assignment G - system optimum and tolls for cost as additional delay
% Compute system optimum with additional delay
clear f;
cvx_begin 
    variable f(M) 
    minimize sum((l.*C).*inv_pos(ones(M,1)-f./C)-l.*C-l.*f)
    subject to 
        B(1,:)*f == 3329
        B(13,:)*f == -3329
        B(2:12,:)*f == 0
        0 <= f <= C
cvx_end

% Social optimum flow with additional delay
f_SocialOptimum_Add = f;

% Compute Wardrop equilibrium with tolls with additional delay
omega = f_SocialOptimum_Add.*(l./C)./((ones(M,1)-f_SocialOptimum_Add./C).^2);
clear f;
cvx_begin
    variable f(M)   
    minimize sum(-l.*C.*log(1-f./C)+f.*omega-l.*f)
    subject to 
        B(1,:)*f == 3329
        B(13,:)*f == -3329
        B(2:12,:)*f == 0
        0 <= f <= C 
cvx_end 

% Wardrop Equilibrium flow with tolls with additional delay
f_Tolls_Add = f;