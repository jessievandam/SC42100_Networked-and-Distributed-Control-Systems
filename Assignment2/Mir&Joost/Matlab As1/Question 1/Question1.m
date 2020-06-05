% Joost Jeschke (4309111) & Miranda van Duijn (4355776)
% Network and Distributed Control SC42100
% Question 1: The rise of the Medici

clear all; close all; clc;

% weight matrix
W = [0 1 1 0 1 0 0 0 0 0 0 0 0 0 0;
     1 0 1 0 0 1 0 0 0 0 0 0 0 0 0
     1 1 0 1 0 1 0 0 0 0 0 0 0 0 0
     0 0 1 0 0 0 0 1 0 0 0 0 0 0 0
     1 0 0 0 0 0 0 1 0 0 0 0 0 0 0
     0 1 1 0 0 0 0 0 0 1 0 0 0 0 0
     0 0 0 0 0 0 0 1 0 1 0 0 0 0 0
     0 0 0 1 1 0 1 0 1 0 1 1 0 0 0
     0 0 0 0 0 0 0 1 0 0 0 0 0 0 0
     0 0 0 0 0 1 1 0 0 0 1 0 1 0 0
     0 0 0 0 0 0 0 1 0 1 0 0 0 1 0
     0 0 0 0 0 0 0 1 0 0 0 0 0 0 1
     0 0 0 0 0 0 0 0 0 1 0 0 0 0 0
     0 0 0 0 0 0 0 0 0 0 1 0 0 0 0
     0 0 0 0 0 0 0 0 0 0 0 1 0 0 0];

% Defining graph
S = [1 1 1 2 2 3 3 4 5 6  7 7  8 8  8  10 10 11 12];
T = [2 3 5 3 6 4 6 8 8 10 8 10 9 11 12 11 13 14 15];
nodenames = {'Castellani' 'Peruzzi' 'Strozzi' 'Ridolfi' 'Barbadori' 'Bischeri' 'Tornabuoni' 'Medici' 'Accaiaiuoli' 'Guadagni' 'Albizzi' 'Salviati' 'Lamberteschi' 'Gimori' 'Pazzi'}';
G = graph(S,T,[],nodenames);
plot(G)

% Graph measures
n = size(W,1); % number of nodes
w = sum(W,2);  % out degree vector
D = diag(w);   % D matrix
P = inv(D)*W;  % normalized weight matrix

%% Exercise 1a: Calculating Bonacich centrality
[V,E] = eig(P');
bonacich = V(:,1)/sum(V(:,1));

%% Exercise 1b: Calculating closeness centrality
dist = distances(G);
closeness = zeros(n,1);
for i = 1:n 
    closeness(i,1) = n/sum(dist(i,:),2);
end
clear i

%% Exercise 1c: Calculating decay centrality
delta = [0.25 0.5 0.75]; %delta values
pi = zeros(n,length(delta)); %initialize pi

for k = 1:length(delta)
    for i = 1:n
        for j = 1:n
            if i == j
            else
                pi(i,k) = pi(i,k) + delta(k)^dist(i,j);
            end
        end
    end
end
clear i j k

%% Exercise 1d: Calculating betweenness centrality
% Calculating all possible paths between all nodes.
% This code is canceled out since it takes quite long to run.
% The paths are stored in AllPaths.mat.
% The code for calculating all the shortest paths from all nodes was found
% on github and copied. The betweenness calculations that follows is ours.

% for i = 1:n
%     for j = 1:n
%         if i ~= j
%             allpaths{i,j} = pathbetweennodes(W,i,j);
%         end
%     end
% end
% clear i j
load AllPaths.mat

%% Calculate the amount of shortest paths per node pair
for i = 1:n
    for j = 1:n
        if i ~= j
            %calculate the length of the shortest paths
            for k = 1:length(allpaths{i,j})
                lengthshortest(k) = length(allpaths{i,j}{k});
            end
            %sort from short to long
            [order, ID] = sort(lengthshortest); 
            %check how many shortest paths there are and save
            for l = 2:length(order)
                if order(l) > order(l-1)
                    ID = ID(1:l-1);
                    order = order(1:l-1);
                    break
                end
            end
            %put shortest paths back in matrix
            for m = 1:length(order)
                allpathsshort{i,j}{m} = allpaths{i,j}{ID(m)};
            end
            clear ID order lengthshortest
        end
    end
end

%% Calculate betweenness
clear i j k l m
% if the length of an element is 2, it cannot contribute to the betweenness
% so erase all the shortest paths of length 2
for i = 1 : n
    for j = 1 : n
        for k = 1 : length(allpathsshort{i,j})
            if length(allpathsshort{i,j}{k}) == 2
                allpathsshort{i,j} = [];
            end
        end
    end
end

%% Betweenness
%betweenenss for node k
bet = zeros(n,1);
for k = 1 : n
    %check over the matrix (i,j)
    for i = 1 : n
        for j = 1 : n
            temp = 0; %to store the times k is in the shortest path
            %check for the amount of elements at (i,j)           
            if length(allpathsshort{i,j}) >= 1
                for m = 1 : length(allpathsshort{i,j})
                    %calculate betweenness for this node combination
                    %exlude starting node and end node 2 : length - 1
                    for l = 2 : length(allpathsshort{i,j}{m})-1 
                        if k == allpathsshort{i,j}{m}(l) %k is in shortest path
                            temp = temp + 1;
                        end
                    end
                end
            %add all the betweenness for a node together
            bet(k) = bet(k) + temp/length(allpathsshort{i,j});
            end
        end
    end
end

%normalize the betweenness values
betn = bet./n^2;

%% Betweenness with the Matlab command for comparison
betweennessvalues = zeros(n,1);
betweennessvalues = centrality(G,'betweenness');
betweennessvaluesnormalized = betweennessvalues/n^2;