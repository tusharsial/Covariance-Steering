close all; clear; clc;

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% Problem data

nx = 1; % dimension of state
I = eye(nx);

% state cost weight matrix
Q = I;

% system matrices
% Example 1
 A = [zeros(nx,1) I(:,1:nx-1)]; 
 B = I(:,end);

%A = rand(nx); B = rand(nx,1);

%% State transition Matrix
M = [A -B*B'; -Q -A'];

% state transition matrix of M
Param.Phi = expm(M);
Param.myu_0 = zeros(nx,1);

% Generate randomly
%Param.Sigma_0 = GenRandomPosDef(nx); % initial covariance
%Param.Sigma_d = GenRandomPosDef(nx); % desired terminal covariance
Param.Sigma_0 = 3;
Param.Sigma_d = 0.5;


%% Fixed point recursion for P0
max_iter = 500; % number of iterations
tol = 1e-8; % tolerance

% initialize
m = 2001; % number of initial guesses
P0 = cell(max_iter,m);
P0_iterates = cell(max_iter,m);

nsym = nx*(nx+1)/2; % Number of distinct elements in P0

% scaling factor for generating random symmetric matrices
a = -5; 
b = 2;

% extract blocks of Phi
Phi = Param.Phi; 
Phi_11 = Phi(1:nx,1:nx); Phi_12 = Phi(1:nx,nx+1:end);
Phi_21 = Phi(nx+1:end,1:nx); Phi_22 = Phi(nx+1:end,nx+1:end);

P = -10:0.01:10;
%P = inv(Param.Sigma_0) - Phi_11/Phi_12;
F = zeros(m,1);

%F = FixedPointRecursionP(Param, A, B, Q, 1.00000000001*P);  

for i=1:m
    F(i) = FixedPointRecursionP(Param, A, B, Q, P(i));   
end  

% plot
figure;
plot(P,F)
hold on;
plot(P,P)
xlabel('P_0')
ylabel('F(P_0)')