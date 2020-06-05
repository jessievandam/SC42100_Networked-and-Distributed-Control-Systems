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
cen_closeness = zeros(n,1); % initialize vector

for i = 1:n 
    cen_closeness(i,1) = n/sum(d(i,:)); % divide number of nodes with sum of distances of row i
end

%% 1C - Decay centrality
pi_d1 = zeros(n,1); % initialize vector
pi_d2 = zeros(n,1); % initialize vector
pi_d3 = zeros(n,1); % initialize vector

for i = 1:n
    for j = 1:n
        if i == j
        else
            % delta = 0.25
            pi_d1(i) = pi_d1(i) + 0.25^dist(i,j);
            
            % delta = 0.5
            pi_d2(i) = pi_d2(i) + 0.5^dist(i,j);
            
            % delta = 0.75
            pi_d3(i) = pi_d3(i) + 0.75^dist(i,j);
        end
    end
end

cen_decay = [pi_d1 pi_d2 pi_d3];

%% 1D - Betweenness centrality
% Compute all posible paths using function found on Github
for i = 1:n
    for j = 1:n
        if i ~= j
            allpaths{i,j} = pathbetweennodes(W,i,j);
        end
    end
end

% Compute distance of shortest path between two nodes
[short_dist] = graphallshortestpaths(sparse(W));

% Compute gij: the fraction of all minimum-distance paths from i to j 
gij = zeros(n,n); % initialization gij
for i = 1:n
    for j = 1:n
        path_ind{i,j} = 0; % initialization path indices per node combination i j
        count = 0;         % initialization counter 1
        countAm = 0;      % initialization counter Amsterdam
        countSch = 0;     % initialization counter Schiphol
        countRot = 0;     % initialization counter Rotterdam
        countUt = 0;      % initialization counter Utrecht
        
        for k = 1:length(allpaths{i,j})
            
            % Compute all path lengths from node i to j
            path_length(k) = length(allpaths{i,j}{k});
            
            % Find the indices of all shortest paths from node i to j
            if (path_length(k)-1)==short_dist(i,j)
                path_ind{i,j} = [path_ind{i,j}; k];
            end
            
            % Check amount of shortest paths from node i to j
            for h = 1:length(path_ind{i,j})
                if (path_ind{i,j}(h,1)>0)
                    count = count+1;  % count 1 if an index is found and there is thus a shortest path
                end
            end
            
            % Check amount of shortest paths that go through a certain node, but not start or end there
            C = cell2mat(allpaths{i,j}(k));
            for h = 1:length(path_ind{i,j})
                
                % Node 4 (Amsterdam)
                if (path_ind{i,j}(h,1)>0) & (ismember(4,C)) & (C(1)~=4) & (C(end)~=4)
                    countAm = countAm+1;
                end
                
                % Node 12 (Schiphol)
                if (path_ind{i,j}(h,1)>0) & (ismember(12,C)) & (C(1)~=12) & (C(end)~=12)
                    countSch = countSch+1;
                end
                
                % Node 11 (Rotterdam)
                if (path_ind{i,j}(h,1)>0) & (ismember(11,C)) & (C(1)~=11) & (C(end)~=11)
                    countRot = countRot+1;
                end
                
                % Node 13 (Utrecht)
                if (path_ind{i,j}(h,1)>0) & (ismember(13,C)) & (C(1)~=13) & (C(end)~=13)
                    countUt = countUt+1;
                end
            end            
        end
        
        GijAm(i,j) = countAm/count;
        GijSch(i,j) = countSch/count;
        GijRot(i,j) = countRot/count;
        GijUt(i,j) = countUt/count;
        
        clear path_length
        clear count
        clear countAm
        clear countSch
        clear countRot
        clear countUt
    end
end

% Summation over all gij's per city
SumGijAm = 0;
SumGijSch = 0;
SumGijRot = 0;
SumGijUt = 0;
for i = 1:n
    for j = 1:n
        if i~=j
            SumGijAm = SumGijAm + GijAm(i,j);
            SumGijSch = SumGijSch + GijSch(i,j);
            SumGijRot = SumGijRot + GijRot(i,j);
            SumGijUt = SumGijUt + GijUt(i,j);
        end
    end
end

% Compute betweenness centrality
betweenAm = SumGijAm*(1/(n^2));
betweenSch = SumGijSch*(1/(n^2));
betweenRot = SumGijRot*(1/(n^2));
betweenUt = SumGijUt*(1/(n^2));

cen_betweenness = [betweenAm; betweenRot; betweenSch; betweenUt];

%% 1D - Function for finding all paths

function pth = pathbetweennodes(adj, src, snk, verbose)
    %PATHBETWEENNODES Return all paths between two nodes of a graph
    %
    % pth = pathbetweennodes(adj, src, snk)
    % pth = pathbetweennodes(adj, src, snk, vflag)
    %
    %
    % This function returns all simple paths (i.e. no cycles) between two nodes
    % in a graph.  Not sure this is the most efficient algorithm, but it seems
    % to work quickly for small graphs, and isn't too terrible for graphs with
    % ~50 nodes.
    %
    % Input variables:
    %
    %   adj:    adjacency matrix
    %
    %   src:    index of starting node
    %
    %   snk:    index of target node
    %
    %   vflag:  logical scalar for verbose mode.  If true, prints paths to
    %           screen as it traverses them (can be useful for larger,
    %           time-consuming graphs). [false]
    %
    % Output variables:
    %
    %   pth:    cell array, with each cell holding the indices of a unique path
    %           of nodes from src to snk.

    % Copyright 2014 Kelly Kearney

    if nargin < 4
        verbose = false;
    end

    n = size(adj,1);

    stack = src;

    stop = false;

    pth = cell(0);
    cycles = cell(0);

    next = cell(n,1);
    for in = 1:n
        next{in} = find(adj(in,:));
    end

    visited = cell(0);

    pred = src;
    while 1

        visited = [visited; sprintf('%d,', stack)];

        [stack, pred] = addnode(stack, next, visited, pred);
        if verbose
            fprintf('%2d ', stack);
            fprintf('\n');
        end

        if isempty(stack)
            break;
        end

        if stack(end) == snk
            pth = [pth; {stack}];
            visited = [visited; sprintf('%d,', stack)];
            stack = popnode(stack);
        elseif length(unique(stack)) < length(stack)
            cycles = [cycles; {stack}];
            visited = [visited; sprintf('%d,', stack)];
            stack = popnode(stack);  
        end
    end
        
    function [stack, pred] = addnode(stack, next, visited, pred)
        newnode = setdiff(next{stack(end)}, pred);
        possible = arrayfun(@(x) sprintf('%d,', [stack x]), newnode, 'uni', 0);

        isnew = ~ismember(possible, visited);

        if any(isnew)
            idx = find(isnew, 1);
            stack = str2num(possible{idx});
            pred = stack(end-1);
        else
            [stack, pred] = popnode(stack);
        end
    end

    function [stack, pred] = popnode(stack)
        stack = stack(1:end-1);
        if length(stack) > 1
            pred = stack(end-1);
        else
            pred = [];
        end
    end
        
 end