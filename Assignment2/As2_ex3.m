%% SC42100 Assignment 2, exercise 3
% Maxime Croft (4390024) and Jessie van Dam (4395832)
clear all; close all; clc;

%% Defining state matrices and initial vector x0 for all four aircrafts
% State vector is (x,y,xdot,ydot)

A1=[1 0 2 0; 0 1 0 2; 0 0 3 0; 0 0 0 3];
B1=[2 0;0 2;3 0;0 3];
x01=[-10;10;-1;1];

A2=[1 0 3 0; 0 1 0 3; 0 0 7 0; 0 0 0 7];
B2=[3 0; 0 3; 7 0; 0 7];
x02=[10;10;1;1];

A3=[1 0 1 0; 0 1 0 1; 0 0 1.1 0; 0 0 0 1.1];
B3=[1 0; 0 1; 1.1 0; 0 1.1];
x03=[10;-10;1;-1];

A4=[1 0 6 0; 0 1 0 6; 0 0 20 0; 0 0 0 20];
B4=[6 0;0 6;20 0; 0 20];
x04=[-10;-10;-1;-1];

% Defining variables
Tfinal = 5; % end time [s]
umax = 100; % upper limit control energy

%% Creating matrices F and H needed for cost function
% Matrix F
F1 = [eye(4); A1; A1^2; A1^3; A1^4; A1^5];
F2 = [eye(4); A2; A2^2; A2^3; A2^4; A2^5];
F3 = [eye(4); A3; A3^2; A3^3; A3^4; A3^5];
F4 = [eye(4); A4; A4^2; A4^3; A4^4; A4^5];

% Matrix H 
n1 = size(B1,1); % number of rows B
n2 = size(B2,2); % number of columns B

G1 = [zeros(n1, 5*n2);
      B1, zeros(n1, 4*n2);
      A1*B1, B1, zeros(n1, 3*n2);
      A1^2*B1, A1*B1, B1, zeros(n1, 2*n2);
      A1^3*B1, A1^2*B1, A1*B1, B1, zeros(n1, n2);
      A1^4*B1, A1^3*B1, A1^2*B1, A1*B1, B1];

G2 = [zeros(n1, 5*n2);
      B2, zeros(n1, 4*n2);
      A2*B2, B2, zeros(n1, 3*n2);
      A2^2*B2, A2*B2, B2, zeros(n1, 2*n2);
      A2^3*B2, A2^2*B2, A2*B2, B2, zeros(n1, n2);
      A2^4*B2, A2^3*B2, A2^2*B2, A2*B2, B2];
  
G3 = [zeros(n1, 5*n2);
      B3, zeros(n1, 4*n2);
      A3*B3, B3, zeros(n1, 3*n2);
      A3^2*B3, A3*B3, B3, zeros(n1, 2*n2);
      A3^3*B3, A3^2*B3, A3*B3, B3, zeros(n1, n2);
      A3^4*B3, A3^3*B3, A3^2*B3, A3*B3, B3];
  
G4 = [zeros(n1, 5*n2);
      B4, zeros(n1, 4*n2);
      A4*B4, B4, zeros(n1, 3*n2);
      A4^2*B4, A4*B4, B4, zeros(n1, 2*n2);
      A4^3*B4, A4^2*B4, A4*B4, B4, zeros(n1, n2);
      A4^4*B4, A4^3*B4, A4^2*B4, A4*B4, B4];

%% Creating matrices H and f for minimization problem quadprog
H1 = 2*(G1'*G1+eye(10));
H2 = 2*(G2'*G2+eye(10));
H3 = 2*(G3'*G3+eye(10));
H4 = 2*(G4'*G4+eye(10));

% Note: for the final equation f still needs to be transposed
f1 = 2*x01'*F1'*G1;
f2 = 2*x02'*F2'*G2;
f3 = 2*x03'*F3'*G3;
f4 = 2*x04'*F4'*G4;

%% Creating constraints for optimization
alpha = 0.01;         % step size
eps = 0.0001;         % stopping criterium
iter = 40000;         % number of iterations
lambda0 = zeros(4,1); % lambda at t=0

% Linear Constraints
Alin = [kron(eye(5), [1 1]); 
        kron(eye(5), [1 -1]); 
        kron(eye(5), [-1 -1]); 
        kron(eye(5), [-1 1])]; 
    
blin = umax/Tfinal * ones(20,1);

% Options quadprog
options = optimoptions('quadprog', 'Display', 'off');

% Initializion iteration loop
stop = sum((x01-x02).^2+(x02-x03).^2+(x03-x04).^2); % stopping criterium
lambda1 = lambda0;                                   
lambda2 = lambda0; 
lambda3 = lambda0; 
lambda4 = lambda0; 

for k = 1:iter
    % Extending f with parts of Lagrangian based on Lagrange multipliers
    f1lambda = f1 + lambda1(:,k)' * G1((Tfinal)*size(A1,1)+1:end,:);
    f2lambda = f2 - lambda1(:,k)' * G2((Tfinal)*size(A1,1)+1:end,:) + lambda2(:,k)' * G2((Tfinal)*size(A1,1)+1:end,:);
    f3lambda = f3 - lambda2(:,k)' * G3((Tfinal)*size(A1,1)+1:end,:) + lambda3(:,k)' * G3((Tfinal)*size(A1,1)+1:end,:);
    f4lambda = f4 - lambda3(:,k)' * G4((Tfinal)*size(A1,1)+1:end,:);

    % Quadratic programming to obtain inputs u
    [u1(:,k), ~] = quadprog(H1, f1lambda', Alin, blin, [], [], [], [], [], options);
    [u2(:,k), ~] = quadprog(H2, f2lambda', Alin, blin, [], [], [], [], [], options);
    [u3(:,k), ~] = quadprog(H3, f3lambda', Alin, blin, [], [], [], [], [], options);
    [u4(:,k), ~] = quadprog(H4, f4lambda', Alin, blin, [], [], [], [], [], options);

    % Compute final state x(5) after each iteration
    x51(:,k) = A1^5*x01 + G1((Tfinal)*size(A1,1)+1:end,:)*u1(:,k); 
    x52(:,k) = A2^5*x02 + G2((Tfinal)*size(A1,1)+1:end,:)*u2(:,k);
    x53(:,k) = A3^5*x03 + G3((Tfinal)*size(A1,1)+1:end,:)*u3(:,k);
    x54(:,k) = A4^5*x04 + G4((Tfinal)*size(A1,1)+1:end,:)*u4(:,k);
    
    % Update Lagrange multipliers
    lambda1(:,k+1) = lambda1(:,k) + alpha*(x51(:,end)-x52(:,end));
    lambda2(:,k+1) = lambda2(:,k) + alpha*(x52(:,end)-x53(:,end));
    lambda3(:,k+1) = lambda3(:,k) + alpha*(x53(:,end)-x54(:,end));
    
    % Check if stopping criterion is satisfied or if iterations are continued
    stop = sum((x51(:,k)-x52(:,k)).^2+(x52(:,k)-x53(:,k)).^2+(x53(:,k)-x54(:,k)).^2);
    if stop < eps
        break;
    end
end

%% Making plots
T = 1:k; % time vector

figure(1);
% Subplot x-coordinate
subplot(221);
hold on; grid on;
plot(T,x51(1,:), 'LineWidth', 1.5);
plot(T,x52(1,:), 'LineWidth', 1.5);
plot(T,x53(1,:), 'LineWidth', 1.5);
plot(T,x54(1,:), 'LineWidth', 1.5);
legend('aircraft 1', 'aircraft 2', 'aircraft 3', 'aircraft 4');
xlabel('iterations'); ylabel('x-coordinate');
% Subplot y-coordinate
subplot(222);
hold on; grid on;
plot(T,x51(2,:), 'LineWidth', 1.5);
plot(T,x52(2,:), 'LineWidth', 1.5);
plot(T,x53(2,:), 'LineWidth', 1.5);
plot(T,x54(2,:), 'LineWidth', 1.5);
legend('aircraft 1', 'aircraft 2', 'aircraft 3', 'aircraft 4');
xlabel('iterations'); ylabel('y-coordinate');
% Subplot xdot
subplot(223);
hold on; grid on;
plot(T,x51(3,:), 'LineWidth', 1.5);
plot(T,x52(3,:), 'LineWidth', 1.5);
plot(T,x53(3,:), 'LineWidth', 1.5);
plot(T,x54(3,:), 'LineWidth', 1.5);
legend('aircraft 1', 'aircraft 2', 'aircraft 3', 'aircraft 4');
xlabel('iterations'); ylabel('xdot');
% Subplot ydot
subplot(224);
hold on; grid on;
plot(T,x51(4,:), 'LineWidth', 1.5);
plot(T,x52(4,:), 'LineWidth', 1.5);
plot(T,x53(4,:), 'LineWidth', 1.5);
plot(T,x54(4,:), 'LineWidth', 1.5);
legend('aircraft 1', 'aircraft 2', 'aircraft 3', 'aircraft 4');
xlabel('iterations'); ylabel('ydot');
suptitle({'Convergence of 4 aircrafts towards common final state x_f,','for states x, y, xdot and ydot'});
