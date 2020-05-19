%% SC42100 Assignment 1, exercise 3
% Maxime Croft (4390024) and Jessie van Dam (4395832)
clear all; close all; clc;

% Load adjacency matrix W Twitter
load -ascii twitter.mat 
W = spconvert(twitter); 

% Adjust W such that it can be used in further calculation
W = [W zeros(size(W,1), size(W,1)-size(W,2))]; % make W square
w = full(sum(W,2));  % out degree vector: column vector containing sum of each row
w(w == 0) = 1;       % add self loop to the node which have no connections

% Measures of Twitter graph
D = sparse(diag(w)); % degree matrix: out degree vector on diagonal
P = inv(D)*W;        % normalized weight/adjacency matrix

%% 3B - PageRank
beta = 0.2;               % PageRank parameter
mu = ones(size(W,1),1);   % intrinsic centrality
iter = 15;                % number of iterations
pr = zeros(size(W,1),10); % initialize matrix

for i = 1:iter
    pr(:,i+1) = pr(:,i) + beta*((1-beta)*(P'))^(i-1)*mu;
end

% find users with highest PageRank
[pr_ranking, ind] = sort(pr(:,end), 'descend');
pr_users = ind(1:5);

%% 3C - Consensus algorithm

% Choosing stubborn nodes
stub1 = 1; % two nodes with high PageRank
stub2 = 2; % two nodes with high PageRank

% Running consensus algorithm
iter = 1000; % number of iterations
clear y;
y = consensus_algorithm(stub1, stub2, P, iter);

% Selecting nodes that get influenced by the stubborn nodes
nod1 = find(y(:,end)==1);

% Creating plot
figure(1);
grid on; 
plot(0:1:iter, y([4561:1:4570],:));
xlabel('Iterations');
ylabel('Opinions');
title({'Simulation of discrete-time consensus algorithm','for stubborn nodes 1 & 2 with values 0 & 1'});
legend('Node 4561', 'Node 4562', 'Node 4563', 'Node 4564', 'Node 4565','Node 4566', 'Node 4567', 'Node 4568', 'Node 4569', 'Node 4570');
axis([0 50 0 1]);

%% 3D - Consensus algorithm - different stubborn nodes

% Choosing stubborn nodes
stub1 = 5232; % two nodes with low PageRank
stub2 = 5233; % two nodes with low PageRank

% Running consensus algorithm
iter = 1000; % number of iterations
clear y;
y = consensus_algorithm(stub1, stub2, P, iter);

% Selecting nodes that get influenced by the stubborn nodes
nod2 = find(y(:,end)==1);
% gives an empty column vector, so no nodes are influenced by the stubborn nodes

% Creating plot 
figure(2);
grid on; 
plot(0:1:iter, y([500:500:5000],:));
xlabel('Iterations');
ylabel('Opinions');
title({'Simulation of discrete-time consensus algorithm','for stubborn nodes 5232 & 5233 with values 0 & 1'});
legend('Node 500', 'Node 1000', 'Node 1500', 'Node 2000', 'Node 2500','Node 3000', 'Node 3500', 'Node 4000', 'Node 4500', 'Node 5000');
axis([0 iter 0 1]);

%% Making function discrete-time consensus algorithm
function y = consensus_algorithm(stub1, stub2, P, iter)
    % Extract B from P
    B = [P(:,stub1) P(:,stub2)]; % B in VxS

    % Q will become matrix P without rows and columns with indices of stubborn nodes
    Q = P; % Q in VxV

    % Remove stubborn nodes from Q and B
    Q(:,[stub1 stub2]) = []; 
    Q([stub1 stub2],:) = []; % Q in RxR
    B([stub1 stub2],:) = []; % B in RxS

    % Prepare for stepwise computation of y
    v = [0; 1];                     % values of stubborn nodes
    y = zeros(length(Q),iter);      % initialization y: creating empty matrix
    y(:,1) = 0.5*ones(length(Q),1); % substituting place 1 in y with y0

    % Discrete-time consensus algorithm 
    for i = 1:iter
        y(:,i+1) = Q*y(:,i) + B*v;    
    end
end
