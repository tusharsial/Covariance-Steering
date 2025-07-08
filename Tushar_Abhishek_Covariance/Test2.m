close all; clear; clc;

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% Problem data

nx = 2; % dimension of state
I = eye(nx);

% state cost weight matrix
Q = I;

% system matrices
% % Example 1
%  A = [zeros(nx,1) I(:,1:nx-1)]; 
%  B = I(:,end);

%Example 3 (stable controllable matrices)
check = 0;
w1 = -2;
w2 = 5;

while check == 0
    A = generate_stable_diagonal_matrix(w1, w2, nx);
    B = rand(nx, 1);
    if (rank(ctrb(A,B))) == 2 % checking controllability
        check = 1;
        break;
    end
end

%% State transition Matrix
M = [A -B*B'; -Q -A'];

% state transition matrix of M
Param.Phi = expm(M);
Param.myu_0 = zeros(nx,1);

% Generate randomly
 Param.Sigma_0 = 3*GenRandomPosDef(nx); % initial covariance
 Param.Sigma_d = GenRandomPosDef(nx); % desired terminal covariance
%Param.Sigma_0 = [0.239085791490040,0.453674299064263;0.453674299064263,1.589545714849335];
%Param.Sigma_d = [19.073123607638990,4.847707958381365;4.847707958381365,15.597445517396638];


%% Fixed point recursion for P0

% extract blocks of Phi
Phi = Param.Phi; 
Phi_11 = Phi(1:nx,1:nx); Phi_12 = Phi(1:nx,nx+1:end);
Phi_21 = Phi(nx+1:end,1:nx); Phi_22 = Phi(nx+1:end,nx+1:end);

a = -4; 
b = -2;
P0_1 = 7*GenRandomSym(a,b,nx);
P0_2 = -12*GenRandomSym(a,b,nx);

%P0_1 = [13.172779435659294,1.106673414098799;1.106673414098799,4.902438711784924];
%P0_2 = [0.395793901951806,0.833683668097575;0.833683668097575,-4.317987637220629];

%P0_1 = FixedPointRecursionP(Param, A, B, Q, P0);
F1 = FixedPointRecursionP(Param, A, B, Q, P0_1);
F2 = FixedPointRecursionP(Param, A, B, Q, P0_2);
% rho = norm(F2 - F1, 2) / norm(P0_2 - P0_1, 2);
disp(eig(F1));
disp(eig(F2));


%P0_1 = FixedPointRecursionP(Param, A, B, Q, P0);

% max_iter = 50;
% log_deltas = zeros(1, max_iter);
% d = zeros(max_iter, 1);
% 
% for k = 2:max_iter
%     P0_2 = FixedPointRecursionP(Param, A, B, Q, P0_1);
%     d(k) = norm(P0_2 - P0_1, 'fro') / norm(P0_1 - P0, 'fro'); % contraction ratio
%     P0 = P0_1;
%     P0_1 = P0_2;
%     disp(P0_1)
% end
% plot(2:max_iter, d(2:end));
% ylabel('Estimated contraction ratio \rho_k');
% xlabel('Iteration');

% log_deltas = zeros(1, 30);
% for k = 1:30
%     P0_2 = FixedPointRecursionP(Param, A, B, Q, P0_1);
%     delta = norm(P0_2 - P0_1, 'fro');
%     log_deltas(k) = log(delta);  % Natural log
%     P0 = P0_1;
%     P0_1 = P0_2;
% end
% 
% plot(1:30, log_deltas, 'o-')
% xlabel('Iteration k')
% ylabel('X_{k+1})')
% title('Convergence of fixed point iteration')
% grid on

