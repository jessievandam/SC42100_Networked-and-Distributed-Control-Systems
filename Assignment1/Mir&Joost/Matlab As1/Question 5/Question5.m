% Joost Jeschke (4309111) & Miranda van Duijn (4355776)
% Network and Distributed Control SC42100
% Question 5: Traffic in LA

clear all; close all;   

% Load data
load -ascii traffic.mat
load -ascii flow.mat
load -ascii traveltime.mat
load -ascii capacities.mat

%% construct graph G
% s and t vectors with the tails and heads of the edges
s = [];
t = [];
% if there is a -1, add to t, if there is a 1, add to s
for i = 1 : size(traffic, 2)
    for j = 1 : size(traffic, 1)
        if traffic(j,i) == 1
            s = [s; j];
        end
        if traffic(j,i) == -1
            t = [t; j];
        end
    end      
end

%% Calculate shortest path between 1 and 17
% create G as a sparse matrix, with the traveltime as weights
G = sparse(s, t, traveltime);
G = [G; zeros(1, size(G,2))]; %make G square

% OUR METHOD
G2 = digraph(s,t);  % defining directed graph G without weights
W2 = adjacency(G2,traveltime);

% show the graph
SP = view(biograph(G,[],'ShowWeights','on'));

% calculate shortest path from 1 to 17
[dist,path,pred]  = graphshortestpath(G,1,13);

% OUR METHOD
[dist2,path2,pred2]  = graphshortestpath(W2,1,13);

%%
% show shortest path by making the nodes and edges red
set(SP.Nodes(path),'Color',[1 0.4 0.4])
edges = getedgesbynodeid(SP,get(SP.Nodes(path),'ID'));
set(edges,'LineColor',[1 0 0])
set(edges,'LineWidth',1.5)

%% Calculate maximum flow between 1 and 17
% create G as a sparse matrix with the flow as weights
Gflow = sparse(s, t, flow);
Gflow = [Gflow; zeros(1, size(Gflow,2))];

% show the graph
MF = view(biograph(Gflow,[],'ShowWeights','on'));

% calculate max flow
[M,F,K] = graphmaxflow(Gflow,1,17);

% show the max flow graph for comparison
MF_2 = view(biograph(F, [], 'ShowWeights', 'on'));

%% Compute inflow and outflow for each node
inout = traffic*flow;

%% Social Equilibrium
lambda = [inout(1); zeros(size(traffic,1)-1,1)];
mu = [zeros(size(traffic,1)-1,1); inout(1)];
dim = size(traffic, 2);
l = traveltime;
B = traffic;
C = capacities;

cvx_begin
    variable f(dim)
    minimize sum((l.*C).*inv_pos(ones(dim,1)-f./C)-l.*C)
    subject to
            B*f == lambda - mu;
            0 <= f <= C
cvx_end

%calculate delay
delay(1) = sum(f.*l.*C./(C-f));

%% Wardrop Equilibrium
cvx_begin
    variable fw(dim)
    minimize sum(-l.*C.*log(C-fw))
    subject to
            B*fw == lambda - mu;
            0 <= fw <= C
cvx_end

%calculate delay
delay(2) = sum(fw.*l.*C./(C-fw));

%% Tolls
cvx_begin
    variable fopt(dim)
    minimize sum(-l.*C.*log(C-fopt)+fopt.*f.*C.*l./((C-f).^2))
    subject to
            B*fopt == lambda - mu;
            0 <= fopt <= C
cvx_end

%calculate delay and weights
delay(3) = sum(fopt.*l.*C./(C-fopt));
tolls(:,1) = f.*C.*l./((C-f).^2);

%% 5g Additional delay
cvx_begin
    variable fg(dim)
    minimize sum((l.*C).*inv_pos(ones(dim,1)-fg./C)-l.*C-l.*fg)

    subject to
            B*fg == lambda - mu;
            0 <= fg <= C
cvx_end

%calculate delay
delay(4) = sum(fg.*l.*C./(C-fg)-l);

%% Additional delay Wardrop
cvx_begin
    variable fgw(dim)
    minimize sum(-l.*C.*log(C-fgw)-l.*fgw)
    subject to
            B*fgw == lambda - mu;
            0 <= fgw <= C
cvx_end

%calculate delay
delay(5) = sum(fgw.*l.*C./(C-fgw)-l);

%% Additional delay Tolls
cvx_begin
    variable fgwt(dim)
    minimize sum(-l.*C.*log(C-fgwt)-l.*fgwt+fgwt.*fg.*C.*l./((C-fg).^2))
    subject to
            B*fgwt == lambda - mu;
            0 <= fgwt <= C
cvx_end

%calculate delay and weights
delay(6) = sum(fgwt.*l.*C./(C-fgwt)-l);
tolls(:,2) = fg.*C.*l./((C-fg).^2);
