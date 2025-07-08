clc;
clear;
close all;

nx = 2;

Sig_d = 13*GenRandomPosDef(nx);
Sig_0 = 2*GenRandomPosDef(nx);

I = eye(nx);
A = rand(nx); B = rand(nx,1);
Q = 0.4*GenRandomPSD(nx); 
%Q = I;

M = [A -B*B'; -Q -A'];

% state transition matrix of M
Phi = expm(M); 
Phi_11 = Phi(1:nx,1:nx); Phi_12 = Phi(1:nx,nx+1:end);
Phi_21 = Phi(nx+1:end,1:nx); Phi_22 = Phi(nx+1:end,nx+1:end);

% % For F3 > 0 
% lamd_max = max(eig(Sig_d));
% c = -5.09/(lamd_max^2);

% % For F3 < 0
% lamd_min = min(eig(Sig_d));
% c = 1.09/(lamd_min^2);

% H1 = zeros(nx);
% 
% P1 = -0.5*(H1 + Sig_d) + sqrtm((0.5*(H1 - Sig_d))^2 + eye(nx));
% F_X = (P1*Phi_12 - Phi_22) \ (Phi_21 - P1*Phi_11);
% 
% H0 = (Phi_21' + Phi_11'*H1)*inv(Phi_22' + Phi_12'*H1);
% 
% X = inv(Sig_0) - H0;


% Test 2
% c = 1;
% X = c*inv(Sig_0);
% X =  -GenRandomPosDef(nx);
X = 

H0 = inv(Sig_0) - X;
H1 = -(Phi_11' - H0*Phi_12') \ (Phi_21' - H0*Phi_22');
P1 = -0.5*(H1 + Sig_d) + sqrtm((0.5*(H1 - Sig_d))^2 + eye(nx));
F_X = (P1*Phi_12 - Phi_22) \ (Phi_21 - P1*Phi_11);

% Display Eigen Values
%disp(eig(X))
disp(eig(X - F_X))