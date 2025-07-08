%% Problem data
close all; clear; clc;

nx = 2; % dimension of state
I = eye(nx);
% state cost weight matrix
Q = I;
% system matrices
%A = [zeros(nx,1) I(:,1:nx-1)]; 
%B = I(:,end);

A = [-1.5 -2; 1 0];
B = [0.5; 0];

disp(rank(ctrb(A,B)));

% block matrix
M = [A -B*B';
    -Q -A'];

% state transition matrix of M
Param.Phi = expm(M);
Param.myu_0 = zeros(nx,1);
%Param.Sigma_0 = GenRandomPosDef(nx); % initial covariance
Param.Sigma_d = GenRandomPosDef(nx); % desired terminal covariance

Param.Sigma_0 =  2*eye(nx);
%Param.Sigma_d = eye(nx)/4; % desired terminal covariance

% unpack the parameters
Phi = Param.Phi; Sig_0 = Param.Sigma_0; Sig_d = Param.Sigma_d;

% number of states
nx = size(Phi,1)/2;

% extract blocks of Phi
Phi_11 = Phi(1:nx,1:nx); Phi_12 = Phi(1:nx,nx+1:end);
Phi_21 = Phi(nx+1:end,1:nx); Phi_22 = Phi(nx+1:end,nx+1:end);

% %% Checking Contractiveness of First Function 
% % initialize
% m = 10; % number of initial guesses
% Y0 = cell(m, 1);
% Sig1 = cell(m, 1);
% frobenius_norms_Y = zeros(nchoosek(m, 2), 1); % Store norms for all pairs of Y0 and Y0 tilde 
% frobenius_norms_Sig = zeros(nchoosek(m, 2), 1); % Store norms for all pairs of Sigma 1 and Sigma 1 tilde
% 
% % Generate Y0 and Sigma 1 values
% for i=1:m
%      Y0{i, 1} = GenRandomSym(nx); % initial guess
%      Sig1{i ,1} = (Phi_12*sqrtm(inv(Sig_0)))*(-eye(nx)/4 + (sqrtm(Sig_0)*( (Sig_0\eye(nx))/2 - Phi_12\Phi_11 - Y0{i,1})*sqrtm(Sig_0) )^2)*(sqrtm(Sig_0\eye(nx))*Phi_12');
% end
% 
% % Compute frobenius norms
% pair_idx = 1;
% for i = 1:m-1
%     for j = i+1:m
%         frobenius_norms_Y(pair_idx, 1) = sqrt(trace((Y0{i, 1} - Y0{j, 1})'*(Y0{i, 1} - Y0{j, 1})));
%         frobenius_norms_Sig(pair_idx, 1) = sqrt(trace((Sig1{i, 1} - Sig1{j, 1})'*(Sig1{i, 1} - Sig1{j, 1})));
%         pair_idx = pair_idx + 1;
%     end
% end
% 
% % Plot the frobenius norms
% figure;
% plot(frobenius_norms_Y);
% hold on;
% plot(frobenius_norms_Sig);
% legend('Y0', 'Sigma1')

% Checking Contractiveness of Third Function
% initialize
m = 10; % number of initial guesses
Sig1 = cell(m, 1);
R = cell(m, 1);
frobenius_norms_Sig = zeros(nchoosek(m, 2), 1); % Store norms for all pairs of Sigma 1 and Sigma 1 tilde
frobenius_norms_R = zeros(nchoosek(m, 2), 1); % Store norms for all pairs of R and R tilde

% Generate Sigma 1 and R values
for i=1:m
     Sig1{i, 1} = GenRandomPosDef(nx); % initial guess
     R{i ,1} = (Sig1{i,1}*Phi_12 - Phi_22)\(Phi_21 - Sig1{i,1}*Phi_11);
end

% Compute frobenius norms
pair_idx = 1;
for i = 1:m-1
    for j = i+1:m
        frobenius_norms_Sig(pair_idx, 1) = sqrt(trace((Sig1{i, 1} - Sig1{j, 1})'*(Sig1{i, 1} - Sig1{j, 1})));
        frobenius_norms_R(pair_idx, 1) = sqrt(trace((R{i, 1} - R{j, 1})'*(R{i, 1} - R{j, 1})));
        pair_idx = pair_idx + 1;
    end
end

% Plot the frobenius norms
figure;
plot(frobenius_norms_Sig);
hold on;
plot(frobenius_norms_R);
legend('Sigma1','R')

% % Fixed point recursion
% num_iter = 100; % number of iterations
% Y0_fpr = cell(num_iter,m);
% Y0_iterates = cell(num_iter,m);
% nsym = nx*(nx+1)/2;
% 
% for i=1:m
%     Y0_fpr{1,i} = Y0{i,1}; % initial guess  
%     idx = triu(true(size(Y0_fpr{1,i}))); % logical index array
%     Y0_iterates{1,i} = Y0_fpr{1,i}(idx);
%     % start recursion
%     for k=1:num_iter-1
%         % update Y0
%         Y0_fpr{k+1,i} = FixedPointRecursion(Param,Y0_fpr{k,i});
%         % extract upper triangular entries in a vector
%         Y0_iterates{k+1,i} = Y0_fpr{k+1,i}(idx);
%     end
% end  
% 
% figure;
% for i=1:nsym
%     subplot(1,nsym,i)
%     for j=1:m
%         x = cellfun(@(v) v(i), Y0_iterates(:,j));
%         semilogy(1:num_iter, x,'-o','LineWidth',2)
%         hold on
%     end
% end
% xlabel('iteration index $k$','FontSize',30)
% ylabel('$Y_0$ entries','FontSize',30)
% 
% %% Plotting Frobenius norms of all combinations (mC2)
% % Preallocate to store Frobenius norms
% frobenius_norms_Y_fpr = zeros(k, nchoosek(m, 2)); % Store norms for all pairs of Y0 and Y0 tilde 
% 
% % Iterate through iterations and calculate Frobenius norm for all pairs
% pair_idx = 1;
% for i = 1:m-1
%     for j = i+1:m
%         for k = 1:num_iter
%             % Calculate Frobenius norm between matrices from initial guesses i and j
%             frobenius_norms_Y_fpr(k, pair_idx) = sqrt(trace((Y0_fpr{k,i} - Y0_fpr{k,j})'*(Y0_fpr{k,i} - Y0_fpr{k,j})));
%         end
%         pair_idx = pair_idx + 1;
%     end
% end
% 
% % Plot Frobenius norm trajectories
% figure;
% hold on;
% for i = 1:size(frobenius_norms_Y_fpr, 2)
%     plot(1:num_iter, frobenius_norms_Y_fpr(:, i));
% end
% hold off;
% title('Frobenius Norm Trajectories for $Y_0$');
% xlabel('Iteration $k$');
% ylabel('Frobenius Norm');