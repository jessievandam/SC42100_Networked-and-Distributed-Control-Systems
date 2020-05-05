%% SC42100 Assignment 1, exercise 1
% Maxime Croft (4390024) and Jessie van Dam (4395832)
clear all; close all; clc;

% Creating graph: specifying graph edges (s,t) in node pairs
s = [1 1 1 2 2 3 3 4 4 4 5 5 6 7 9 9 11 11 12];
t = [4 7 8 4 10 4 13 7 8 12 9 11 11 9 12 13 12 13 13];
nodenames = {'Alkmaar','Almere','Amersfoort','Amsterdam','Den Haag','Dordrecht','Haarlem','Hoorn','Leiden','Lelystad','Rotterdam','Schiphol','Utrecht'};
G = graph(s,t,[],nodenames);  % defining graph G
plot(G);                      % plotting graph G

% Weight/adjacency matrix W
W = [0 0 0 1 0 0 1 1 0 0 0 0 0;
    0 0 0 1 0 0 0 0 0 1 0 0 0;
    0 0 0 1 0 0 0 0 0 0 0 0 1;
    1 1 1 0 0 0 1 1 0 0 0 1 0;
    0 0 0 0 0 0 0 0 1 0 1 0 0;
    0 0 0 0 0 0 0 0 0 0 1 0 0;
    1 0 0 1 0 0 0 0 1 0 0 0 0;
    1 0 0 1 0 0 0 0 0 0 0 0 0;
    0 0 0 0 1 0 1 0 0 0 0 1 1;
    0 1 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 1 1 0 0 0 0 0 1 1;
    0 0 0 1 0 0 0 0 1 0 1 0 1;
    0 0 1 0 0 0 0 0 1 0 1 1 0];

% Measures of graph G
n = size(W,1);    % total number of nodes
w = sum(W,2);     % out degree vector: column vector containing sum of each row
D = diag(w);      % degree matrix: out degree vector on diagonal
P = inv(D)*W;     % normalized weight/adjacency matrix
d = distances(G); % matrix d where d(i,j) is the length of the shortest path between node i and node j 

%% 1A - Bonacich centrality
[V,D] = eig(P');              % matrix V columns contain right eigenvectors of P'
pii = V(:,1);                 % eigenvector pi corresponding to an eigenvalue of 1
cen_bonacich = pii/sum(pii);  % normalizing pi to get Bonacich centrality

%% 1B - Closeness centrality
cen_closeness = zeros(n,1);
for i = 1:n 
    cen_closeness(i,1) = n/sum(d(i,:)); % divide number of nodes with sum of distances of row i
end

%% 1C - Decay centrality

%% 1D - Betweenness centrality


