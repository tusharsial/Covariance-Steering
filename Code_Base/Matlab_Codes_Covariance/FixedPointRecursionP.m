%%% This function solves the fixed point recursion algorithm to compute 
% initial Costate matrix P0 for our Optimal Covariance Steering Algorithm

function P0 = FixedPointRecursionP(nx, Param, A, B, Q, P0)

% unpack the system parameters
Sig_0 = Param.Sigma_0; 
Sig_d = Param.Sigma_d;

% extract blocks of Phi
Phi = Param.Phi; 
Phi_11 = Phi(1:nx,1:nx); Phi_12 = Phi(1:nx,nx+1:end);
Phi_21 = Phi(nx+1:end,1:nx); Phi_22 = Phi(nx+1:end,nx+1:end);

%% Recursion
% Step 1: P0 -> H0 (Linear Map) 
H0 = inv(Sig_0) - P0;

% Step 2: H0 -> H1 (Linear Fractional Transformation (LFT))
H1 = -(Phi_11' - H0*Phi_12') \ (Phi_21' - H0*Phi_22');

% Step 3: H1 -> P1 (Solve the algebraic Riccatti Equation)
P1 = -0.5*(H1 + Sig_d) + sqrtm((0.5*(H1 - Sig_d))^2 + eye(nx));

% Step 4: P1 -> P0 (Linear Fractional Transformation (LFT)) 
P0 = (P1*Phi_12 - Phi_22) \ (Phi_21 - P1*Phi_11);

%% (Optional: To check whether the above LFTs are consistent with our Ricatti 
% differential equation in H, I've implemented the integrator for 
% Matrix Ricatti Equation H and P. So you can replace the respective LFT 
% steps with the differential equations given below

%N = 500; % number of grid points

%%% Integrator for Matrix Ricatti Equation H (Forward Integration)
% t_fineH = linspace(t_initial, t_final, N);
% [t,H] = ode45(@(t,H) LyapODE_H(t,H,A,B,Q), t_fineH, H0(:));
% H1 = reshape(H(end,:), size(A));

%%% Integrator for Matrix Ricatti Equation P (Backward Integration)
%t_fineP = linspace(t_final, t_initial, N);
% [t,P] = ode45(@(t,P) LyapODE_P(t,P,A,B,Q), t_fineP, P1(:));
% P0 = reshape(P(end,:), size(A));

%% Additional Checks that you can do 
% You can check the eigenvalues of the matrices at each step to see whether
% our converged solution was indeed feasible without any violations