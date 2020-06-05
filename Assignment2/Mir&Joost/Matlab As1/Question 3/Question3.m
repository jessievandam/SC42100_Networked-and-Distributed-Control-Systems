% Joost Jeschke (4309111) & Miranda van Duijn (4355776)
% Network and Distributed Control SC42100
% Question 3: Influence on Twitter

clear all;

load -ascii users.mat
load -ascii twitter.mat
W = spconvert(twitter);
W = [W zeros(size(W,1),size(W,1)-size(W,2))]; %make W square

w = W*ones(size(W,1),1); %calculate W
w(w == 0) = 1; %add self loop to the node which have no connections
D = sparse(diag(w)); %diagonalize w
P = D \ W; %calculate P

%% Pagerank
%page rank variables
beta = 0.15;
mu = ones(size(W,1),1);

%intialize PageRank
PR = zeros(size(W,1),10);

%iteratively PageRank calculation
for i = 1 : 10
    PR(:,i+1) = PR(:,i) + beta*(1-beta)^(i-1)*(P')^(i-1)*mu;
end

%% Calculate nodes with highest PageRank
[order, ID] = sort(PR(:,end), 'descend');
ID(1:5);

%% Compare with Matlab function
G = digraph(W);
pr = centrality(G, 'pagerank', 'FollowProbability', 0.85);
[order2, ID2] = sort(pr, 'descend');
ID2(1:5);

%% Stubborn nodes
stub = [15 16]; %choose stubborn nodes
u = [1; 0]; %choose which one has value 0 and 1
iter = 1000; %number of iterations
y = concensus(P, stub, u, iter); %concensus function

%plot y (selection of nodes that get influenced by both stubborn nodes)
figure; plot(0:1:iter, y(255:264, :));
xlabel('Iterations');
ylabel('Opinions');
title('Stubborn nodes 15 & 16');
legend('Node 255', 'Node 256', 'Node 257', 'Node 258', 'Node 259', ...
     'Node 260', 'Node 261', 'Node 262', 'Node 263', 'Node 264');
axis([0 iter/10 0 1]);


%% Different stubborn nodes
stub = [1 16];
y1 = concensus(P, stub, u, iter);
figure;
plot(0:1:iter, y1(255:264,:));
xlabel('Iterations');
ylabel('Opinions');
title('Stubborn nodes 1 & 16');
legend('Node 255', 'Node 256', 'Node 257', 'Node 258', 'Node 259', ...
     'Node 260', 'Node 261', 'Node 262', 'Node 263', 'Node 264');
axis([0 iter 0 1]);