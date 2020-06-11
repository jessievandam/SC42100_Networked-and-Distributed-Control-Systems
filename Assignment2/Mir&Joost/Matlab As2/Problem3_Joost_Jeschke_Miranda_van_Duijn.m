%% Joost Jeschke (4309111) & Miranda van Duijn (4355776)
clear all;
close all;
clc;

%% state vector is (x,y,xdot,ydot)
A(:,:,1)=[1 0 2 0; 0 1 0 2; 0 0 3 0; 0 0 0 3];
B(:,:,1)=[2 0;0 2;3 0;0 3];
x0(:,1)=[-10;10;-1;1];

A(:,:,2)=[1 0 3 0; 0 1 0 3; 0 0 7 0; 0 0 0 7];
B(:,:,2)=[3 0; 0 3; 7 0; 0 7];
x0(:,2)=[10;10;1;1];

A(:,:,3)=[1 0 1 0; 0 1 0 1; 0 0 1.1 0; 0 0 0 1.1];
B(:,:,3)=[1 0; 0 1; 1.1 0; 0 1.1];
x0(:,3)=[10;-10;1;-1];

A(:,:,4)=[1 0 6 0; 0 1 0 6; 0 0 20 0; 0 0 0 20];
B(:,:,4)=[6 0;0 6;20 0; 0 20];
x0(:,4)=[-10;-10;-1;-1];

Tfinal=5;
umax=100;

%% Construct F and phi matrices. H and h for optimization
for i = 1 : size(A,3)
    F(:,:,i) = zeros((Tfinal+1)*size(A,1), size(A,2));
    phi(:,:,i) = zeros((Tfinal+1)*size(A(:,:,1)*B(:,:,1),1),(Tfinal)*size(A(:,:,1)*B(:,:,1),2));     
    F(1:size(A,1), 1:size(A,2), i) = eye(4); %first part is eye(4)
    
    for j = 2 : Tfinal+1 %Column
        for k = 1 : Tfinal %Row
            if k < j %B on the diagonal
                phi((j-1)*size(A,1)+1:j*size(A,1), (k-1)*size(B,2)+1:k*size(B,2), i) = B(:,:,i);
            end
            if k < j - 1 %A^(j-k)*B below the diagonal
                phi((j-1)*size(A,1)+1:j*size(A,1), (k-1)*size(B,2)+1:k*size(B,2), i) = A(:,:,i)^(j-k-1)*phi((j-1)*size(A,1)+1:j*size(A,1), (k-1)*size(B,2)+1:k*size(B,2), i);
            end            
        end
    end
    
    for j = 2 : Tfinal+1
        F((j-1)*size(A,1)+1:j*size(A,1), 1:size(A,2), i) = A(:,:,i)^(j-1); %Powers of A below eachother
    end
                 
    
    lagr(:,:,i) = phi((Tfinal)*size(A,1)+1:end, :, i); %The part of Phi for the Lagrangian
    
    % functions for quadprog
    H(:,:,i) = phi(:,:,i)'*phi(:,:,i)+eye(10);
    h(:,:,i) = 2*x0(:,i)'*F(:,:,i)'*phi(:,:,i);
    
    final(:,i) = x0(:,i); %initialize final values for stop
end


%% Optimization
for i = 1 : 3 
    lambda(:,1,i) = 0*ones(4,1); % initial values lambda
end

alpha = 1e-2; % step size

% Linear Constraints
Abound = [kron(eye(5), [1 1]); kron(eye(5), [1 -1]); kron(eye(5), [-1 -1]); kron(eye(5), [-1 1])]; 
Bbound = umax/Tfinal*ones(20,1);

% Do not display for Quadprog
options = optimoptions('quadprog', 'Display', 'off');

stop = sum((final(:,1)-final(:,2)).^2+(final(:,2)-final(:,3)).^2+(final(:,3)-final(:,4)).^2);

for k = 1 : 30000 %iters
    for j = 1 : size(A,3) %number of planes
        %Quadprog 0.5*u'Hu + h'u, add right lambda parts to h
    	hlambda(:,:,1) = h(:,:,1) + lambda(:,k,1)'*lagr(:,:,1);
        hlambda(:,:,2) = h(:,:,2) - lambda(:,k,1)'*lagr(:,:,2) + lambda(:,k,2)'*lagr(:,:,2);
        hlambda(:,:,3) = h(:,:,3) - lambda(:,k,2)'*lagr(:,:,3) + lambda(:,k,3)'*lagr(:,:,3);
        hlambda(:,:,4) = h(:,:,4) - lambda(:,k,3)'*lagr(:,:,4);
        
        %Quadprog, calculate perfect inputs
        [X(:,k,j), Fval(k,j)] = quadprog(2*H(:,:,j), hlambda(:,:,j)', Abound, Bbound, [], [], [], [], [], options);
        
        %Calculate final location
        location(:, k, j) = A(:,:,j)^5*x0(:,j)+lagr(:,:,j)*X(:,k,j); %A^5*x0+A^4B*u(0)+A^3Bu(1) etc.
        final(:,j) = location(:, k, j); %update final values for stop
    end
    
    % Update lambda's
    lambda(:,k+1,1) = lambda(:,k,1)+alpha*(location(:,end,1)-location(:,end,2));
    lambda(:,k+1,2) = lambda(:,k,2)+alpha*(location(:,end,2)-location(:,end,3));
    lambda(:,k+1,3) = lambda(:,k,3)+alpha*(location(:,end,3)-location(:,end,4));
    
    % Display k for to follow iterations
    stop = sum((final(:,1)-final(:,2)).^2+(final(:,2)-final(:,3)).^2+(final(:,3)-final(:,4)).^2);
    
    if stop < 0.0001
        break;
    end
end


%% Plot result
T = 1 : k; % Time vector for plots
figure;
for i = 1 : size(A,3) % Plot state 5 of x, y, dx, dy in 4 different plots
    subplot(2,2,1); hold on; plot(T, location(1,:,i));
    subplot(2,2,2); hold on; plot(T, location(2,:,i));
    subplot(2,2,3); hold on; plot(T, location(3,:,i));
    subplot(2,2,4); hold on; plot(T, location(4,:,i));
end
subplot(2,2,1); legend('Plane 1', 'Plane 2', 'Plane 3', 'Plane 4'); title('x'); xlabel('Iterations'); ylabel('x');
subplot(2,2,2); legend('Plane 1', 'Plane 2', 'Plane 3', 'Plane 4'); title('y'); xlabel('Iterations'); ylabel('y');
subplot(2,2,3); legend('Plane 1', 'Plane 2', 'Plane 3', 'Plane 4'); title('xdot'); xlabel('Iterations'); ylabel('dx');
subplot(2,2,4); legend('Plane 1', 'Plane 2', 'Plane 3', 'Plane 4'); title('ydot'); xlabel('Iterations'); ylabel('dy');

