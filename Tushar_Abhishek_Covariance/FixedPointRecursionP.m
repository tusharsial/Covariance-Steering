function P0 = FixedPointRecursionP(Param, A, B, Q, P0)
%disp(eig(P0));
% unpack the parameters
Sig_0 = Param.Sigma_0; Sig_d = Param.Sigma_d;

nx = 2;

% extract blocks of Phi
Phi = Param.Phi; 
Phi_11 = Phi(1:nx,1:nx); Phi_12 = Phi(1:nx,nx+1:end);
Phi_21 = Phi(nx+1:end,1:nx); Phi_22 = Phi(nx+1:end,nx+1:end);

%% Recursion
t_initial = 0; 
t_final = 1;
N = 500; % number of grid points

H0 = inv(Sig_0) - P0;

%disp(eig(H0));

%% integrator for Matrix Ricatti Equation H
% t_fineH = linspace(t_initial, t_final, N);
% [t,H] = ode45(@(t,H) LyapODE_H(t,H,A,B,Q), t_fineH, H0(:));
% H1 = reshape(H(end,:), size(A));
%P0 = reshape(H(end,:), size(A));
H1 = -(Phi_11' - H0*Phi_12') \ (Phi_21' - H0*Phi_22');

%disp(eig(H1));

%% (Solve the algebraic Riccatti Equation)
%P1 = ARE(H1, Sig_d);
P1 = -0.5*(H1 + Sig_d) + sqrtm((0.5*(H1 - Sig_d))^2 + eye(nx));
%P0 = ARE(P0, Sig_d);
%disp(eig(P1));

%% integrator for Matrix Ricatti Equation P
%N = 500; % number of grid points
%t_fineP = linspace(t_final, t_initial, N);

% [t,P] = ode45(@(t,P) LyapODE_P(t,P,A,B,Q), t_fineP, P1(:));
% P0 = reshape(P(end,:), size(A));
P0 = (P1*Phi_12 - Phi_22) \ (Phi_21 - P1*Phi_11);

% % Test
% N = 500; % number of grid points
% t_fineP = linspace(t_initial, t_final, N);
% 
% [t,P] = ode45(@(t,P) LyapODE_P(t,P,A,B,Q), t_fineP, P0(:));
% P1_test = reshape(P(end,:), size(A));

%disp(H0);
%disp(H1);
%disp(P1);


