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

%% 3C - consensus algorithm

stub = [15 16]; %choose stubborn nodes
u = [1; 0]; %choose which one has value 0 and 1
iter = 1000; %number of iterations

%Extract B from W
    B = [P(:,stub(1)) P(:,stub(2))];
    
%Remove stubborn nodes from W and B
for i = 1 : length(stub)
    P(stub(i),:) = [];
    P(:,stub(i)) = [];
    B(stub(i),:) = [];
end
Q = P; %store in Q

%initialize y
y = zeros(size(Q,1),iter);
y(:,1) = 0.5*ones(size(Q,1),1);

for i = 1 : iter
    y(:,i+1) = Q*y(:,i) + B*u;    
end


%plot y (selection of nodes that get influenced by both stubborn nodes)
figure; plot(0:1:iter, y(255:264, :));
xlabel('Iterations');
ylabel('Opinions');
title('Stubborn nodes 15 & 16');
legend('Node 255', 'Node 256', 'Node 257', 'Node 258', 'Node 259', ...
     'Node 260', 'Node 261', 'Node 262', 'Node 263', 'Node 264');
axis([0 iter/10 0 1]);


