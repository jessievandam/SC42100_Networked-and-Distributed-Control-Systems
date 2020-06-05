clear all
close all
load -ascii twitter.mat
W = spconvert(twitter);
[q,q1] = size(W);
W(:,q1+1:q) = zeros(q,q-q1);

%%
w_in = W*ones(q,1);
w_out = W'*ones(q,1);

%%
mu = ones(q,1);
beta = 0.15;
W1 = W*eye(q);
w = W1*ones(q,1); 
w(w==0) = 1;
D = diag(w);
P = D\W;
P2 = P';
%%
n = 10 ;
PageRank = zeros(q,1);
for i = 1:n
    PageRank = PageRank + beta*(((1-beta)^(i-1))*(P2)^(i-1))*mu;
end
%%
descending_order = sort(PageRank(:,1),'descend');

for i = 1:5
    central_nodes(i) = find(PageRank(:,1) == descending_order(i),1);
end

%%
stubborn1 = 11;
stubborn2 = 12;
iterations = 1000;
clear y
y = Stubborn_nodes(stubborn1,stubborn2,P,iterations);
figure(1)
plot(y(1:10,:)')
legend('node1','node2','node3','node4','node5','node6','node7','node8','node9','node10')
axis([1, iterations, 0, 1])


%%
stubborn1 = 1;
stubborn2 = 2;
iterations = 1000;
clear y
y = Stubborn_nodes(stubborn1,stubborn2,P,iterations);
figure(2)
plot(y')
axis([1, iterations, 0, 1])



%%
function y = Stubborn_nodes(stubborn1,stubborn2,P,iterations)

Q = P;
Q([stubborn1 stubborn2],:) = [];

B = [ Q(:,stubborn1), Q(:,stubborn2) ];
u = [ 0 ; 1 ];

Q(:,[stubborn1 stubborn2]) = [];

y = zeros(length(Q(:,1)),iterations);
y = ones(length(Q(:,1)),1);

for i = 1:iterations
    y(:,i+1) = Q*y(:,i) + B*u;
end
end

