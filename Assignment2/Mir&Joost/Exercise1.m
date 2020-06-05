%% SC4210 NDC: Assignment 1, Exercise 1
% Joost Jeschke (4309111) and Miranda van Duijn (4355776)
clear all; close all; clc;

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
     0 0 0 0 0 0 0 0 0 0 0 1 0 0 0]; % weight matrix

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

% Graph measures via Matlab
Wmat = adjacency(G);
wmat = Wmat*ones(height(G.Nodes),1);
Dmat = diag(wmat);
Pmat = inv(Dmat)*Wmat;

%% Exercise 1a: Calculating Bonacich centrality
[V,E] = eig(P');
bonacich = V(:,1)/sum(V(:,1));

[Vmat,Emat] = eig(Pmat');
bonacichmat = Vmat(:,1)/sum(Vmat(:,1));

%% Exercise 1b: Calculating closeness centrality
dist = distances(G);
closeness = zeros(n,1);
for i = 1:n 
    closeness(i,1) = n/sum(dist(i,:),2);
end
clear i

%% Exercise 1c: Calculating decay centrality
delta = [0.25 0.5 0.75];
pi = zeros(n,length(delta));

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
betweennessvalues = zeros(n,1);

% Find all shortest paths
for i = 1:n
    [R C] = find(dist == 1)

end

betweennessvalues = centrality(G,'betweenness');
betweennessvaluesnormalized = betweennessvalues/(n^2);

%% Calculating all possible paths between all nodes
for i = 1:n
    for j = 1:n
        if i ~= j
            allpaths{i,j} = pathbetweennodes(W,i,j);
        end
    end
end
clear i j

%% 
for i = 1:n
    allpaths{i,i} = []; % paths were returned from node to same node so clear these cells
end
clear i j

for i = 1:n
    for j = 1:n
        for k = 1:length(allpaths{i,j})
            lengthshortest = allpaths{i,j}
%             shortestpaths{i,j} = ;
        end
    end
end


