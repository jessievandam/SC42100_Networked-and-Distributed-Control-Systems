clear all
close all

n = 15;
W = zeros(n,n); %Creating empty W matrix

%Adding the relevant families
Medici = 12;
Strozzi = 6;
Tornabuoni = 7;

%Adding weights from the links between families
W(1,2) = 1;
W(1,6) = 1;
W(1,11) = 1;

W(2,3) = 1;
W(2,6) = 1;


W(3,4) = 1;
W(3,6) = 1;

W(4,5) = 1;
W(4,7) = 1;
W(4,9) = 1;

W(6,8) = 1;

W(7,12) = 1;

W(8,12) = 1;

W(9,10) = 1;
W(9,12) = 1;

W(11,12) = 1;

W(12,13) = 1;
W(12,14) = 1;

W(14,15) = 1;

W = W + W'; %Summing with its transposed as the graph is undirected, with all weights on the links equal to 1

G = digraph(W);
figure
plot(G)
%% Finding the Bonacich centrality vector
w = W*ones(n,1);

D = diag(w);
P = D\W;

[eigvector,eigvalue] = eig(P');
%% Taking the Bonacich centrality for the different famlilies Medici, Strozzi an Tornabuoni.
Bona = [eigvector(Medici,1), eigvector(Strozzi,1), eigvector(Tornabuoni,1)];

%% Calculating the closeness centrality for the different famlilies Medici, Strozzi an Tornabuoni.
closeness = zeros(1,3);

%Summing the shortest path to every different node
for i = 1:n
    closeness(1,1) = length(shortestpath(G,Medici,i))-1 + closeness(1,1);
    closeness(1,2) = length(shortestpath(G,Strozzi,i))-1 + closeness(1,2);
    closeness(1,3) = length(shortestpath(G,Tornabuoni,i))-1 + closeness(1,3);
end


%closeness centrality
closeness = 1./(closeness/n);

%% Calculating the decay centrality for the different famlilies Medici, Strozzi an Tornabuoni.

delta = 0.25;

%Creating the decay centrality vectors for 3 different delta's.
%Pre-sustracting 1 to remove the effect from the distance to itself.
decay1 = -ones(3,1);
decay2 = -ones(3,1);
decay3 = -ones(3,1);

%Summing the decay centrality to all nodes (including itself).
for i = 1:n
    for j = 1:3
        decay1(j,1) = decay1(j,1) + (j*delta)^(length(shortestpath(G,Medici,i))-1);
        decay2(j,1) = decay2(j,1) + (j*delta)^(length(shortestpath(G,Strozzi,i))-1);
        decay3(j,1) = decay3(j,1) + (j*delta)^(length(shortestpath(G,Tornabuoni,i))-1);
    end
end


%%
k = zeros(n,1);
for i = 1:n
    for j = 1:n
        %         if j>i
        G1 = G;
        Paths = shortestpath(G1,i,j);
        spath = Paths;
        L = length(Paths)-1;
        
        for p = 2:L
            G1 = G;
            G1 = rmnode(G1,Paths(p));
            
            ii = i;
            iii = j;
            if i > Paths(p)
                ii = ii-1;
            end
            if iii > Paths(p)
                iii = iii -1;
            end
            L1 = length(shortestpath(G1,iii,ii))-1;
            if L1 == L
                spath(p,:) = shortestpath(G1,iii,ii);
                
                
                Medici = 12;
                Strozzi = 6;
                Tornabuoni = 7;
                if Medici > Paths(p)
                    Medici = Medici - 1;
                end
                if Strozzi > Paths(p)
                    Strozzi = Strozzi - 1;
                end
                if Tornabuoni > Paths(p)
                    Tornabuoni = Tornabuoni - 1;
                end
                k1(j,i) = sum(spath(:) == Medici);
                k2(j,i) = sum(spath(:) == Strozzi);
                k3(j,i) = sum(spath(:) == Tornabuoni);
            end
        end
        
    end
    %     end
end



%%
minpaths = graphallshortestpaths(sparse(W));
k0 = zeros(n,n);
k1 = zeros(n,n);
k2 = zeros(n,n);
k3 = zeros(n,n);

for i = 1:15
    for j = 1:15
        if i ~= j
            Paths = All_Paths(W, i, j, n);
            for k = 1:length(Paths)
                pathtaken = cell2mat(Paths(k));
                if length(pathtaken) -1 == minpaths(i,j)
                    k0(i,j) = k0(i,j) + 1;
                    if any(pathtaken == Medici)
                        k1(i,j) = k1(i,j) + 1;
                    end
                    if any(pathtaken == Strozzi)
                        k2(i,j) = k2(i,j) + 1;
                    end
                    if any(pathtaken == Tornabuoni)
                        k3(i,j) = k3(i,j) + 1;
                    end
                end
            end
        end
    end
end

k0(k0 == 0) = Inf;
betweenness(1) = (1/(n)^2)*(sum(sum(k1./k0))+1);
betweenness(2) = (1/(n)^2)*(sum(sum(k2./k0))+1);
betweenness(3) = (1/(n)^2)*(sum(sum(k3./k0))+1);

function Paths = All_Paths(W, Source, Sink, n)
%Create a list of all connected nodes
for in = 1:n
    Connected_Nodes{in} = find(W(in,:));
end
Stack_Nodes = Source;
Next_Node = Source;
Paths = cell(0);
Cycle = cell(0);
Visited = cell(0);

while 1
    %Update the list of visited nodes
    Visited = [Visited; sprintf('%d,', Stack_Nodes)];
    [Stack_Nodes, Next_Node] = Add_Node(Stack_Nodes, Connected_Nodes, Visited, Next_Node);
    
    if isempty(Stack_Nodes)
        break;
    end
    
    if Stack_Nodes(end) == Sink
        Paths = [Paths; {Stack_Nodes}];
        Visited = [Visited; sprintf('%d,', Stack_Nodes)];
        Stack_Nodes = Remove_Node(Stack_Nodes);
    elseif length(unique(Stack_Nodes)) < length(Stack_Nodes)
        Cycle = [Cycle; {Stack_Nodes}];
        Visited = [Visited; sprintf('%d,', Stack_Nodes)];
        Stack_Nodes = Remove_Node(Stack_Nodes);
    end
end
    
    function [Stack_Nodes, Next_Node] = Add_Node(Stack_Nodes, Connected_Nodes, Visited, Next_Node)
        New_Nodes = setdiff(Connected_Nodes{Stack_Nodes(end)}, Next_Node);                    %Sort all nodes that are not the next node
        Possible_Node = arrayfun(@(x) sprintf('%d,', [Stack_Nodes x]), New_Nodes, 'uni', 0);  %Specity which nodes are possible
        New_Node = ~ismember(Possible_Node, Visited);                                   %
        
        if any(New_Node)
            Stack_Nodes = str2num(Possible_Node{find(New_Node,1)});
            Next_Node = Stack_Nodes(end-1);
        else
            [Stack_Nodes, Next_Node] = Remove_Node(Stack_Nodes);
        end
    end
    
    function [Stack_Nodes, Next_Node] = Remove_Node(Stack_Nodes)
        Stack_Nodes = Stack_Nodes(1:end-1);                 %Remove last node from the stack
        if length(Stack_Nodes) > 1
            Next_Node = Stack_Nodes(end-1);           
        else
            Next_Node = [];                     %Last node has been reached
        end
    end
end