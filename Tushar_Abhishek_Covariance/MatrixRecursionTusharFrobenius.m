close all; clear; clc;

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% Problem data

nx = 2; % dimension of state

% scaling factor for generating random symmetric matrices
a = -5; 
b = 2;

I = eye(nx);

% state cost weight matrix
Q = GenRandomPosDef(nx);
%Q = I;

% system matrices
% % Example 1
A = [zeros(nx,1) I(:,1:nx-1)]; 
B = I(:,end);

%A = rand(nx); B = rand(nx,1);

% % Example 2
% A = [-1.5 -2; 7 0];
% B = [0.5; 0];

% %Example 3 (stable controllable matrices)
% check = 0;
% w1 = -2;
% w2 = 5;
% 
% while check == 0
%     A = generate_stable_diagonal_matrix(w1, w2, nx);
%     B = rand(nx, 1);
%     if (rank(ctrb(A,B))) == 2 % checking controllability
%         check = 1;
%         break;
%     end
% end
 
%% %Example 4 (Linearized CW equations)
% % The satellite is moving in a LEO (Low earth Orbit) corresponding to an 
% % orbital period of about 93 minutes;
% a = 6793237; % semimajor axis (in m)
% mu_e = 3.986*10^14; % in m^3/s^2
% n = sqrt(mu_e/a^3); % orbital rate (in m/s)
% 
% % Define subblocks for the CW equations
% A11 = zeros(nx/2);
% A12 = eye(nx/2);
% A21 = [3*n^2 0 0; 0 0 0; 0 0 -n^2];
% A22 = [0 2*n 0; -2*n 0 0; 0 0 0];
% 
% A = [A11 A12; A21 A22]; 
% B = [A11; A12];

%% State transition Matrix
M = [A -B*B'; -Q -A'];

% state transition matrix of M
Param.Phi = expm(M);
Param.myu_0 = zeros(nx,1);

% Generate randomly
Param.Sigma_0 = 5*GenRandomPosDef(nx); % initial covariance
Param.Sigma_d = GenRandomPosDef(nx); % desired terminal covariance
%Param.Sigma_0 = [0.239085791490040,0.453674299064263;0.453674299064263,1.589545714849335];
%Param.Sigma_d = [19.073123607638990,4.847707958381365;4.847707958381365,15.597445517396638];

% % Specific Example
% Param.Sigma_0 =  2*eye(nx); % initial covariance
% Param.Sigma_d = eye(nx); % desired terminal covariance

% % Displaying initial distance 
% disp(sqrt(trace((Param.Sigma_0 - Param.Sigma_d)'*(Param.Sigma_0 - Param.Sigma_d)))/2);

%% Fixed point recursion for P0
max_iter = 500; % number of iterations
tol = 1e-8; % tolerance

% initialize
m = 1; % number of initial guesses
P0 = cell(max_iter,m);
P0_iterates = cell(max_iter,m);

% % Fill every cell with a zero matrix of size nx × nx
% for i = 1:max_iter
%     for j = 1:m
%         P0{i,j} = zeros(nx); % nx × nx zero matrix
%     end
% end

nsym = nx*(nx+1)/2; % Number of distinct elements in P0

% extract blocks of Phi
Phi = Param.Phi; 
Phi_11 = Phi(1:nx,1:nx); Phi_12 = Phi(1:nx,nx+1:end);
Phi_21 = Phi(nx+1:end,1:nx); Phi_22 = Phi(nx+1:end,nx+1:end);

for i=1:m
    P0{1,i} = GenRandomSym(a,b,nx); % initial guess
    %P0{1,i} = inv(Param.Sigma_0) - Phi_11/Phi_12; % initial guess
    %P0{1,i} = [13.172779435659294,1.106673414098799;1.106673414098799,4.902438711784924]; % initial guess
    idx = triu(true(size(P0{1,i}))); % logical index array
    P0_iterates{1,i} = P0{1,i}(idx);
    
    % Declaring error and iterate index
    err = 1000;
    r = 1;

    % start recursion
    while err > tol && r<max_iter-1
        % update Y0
        P0{r+1,i} = FixedPointRecursionP(Param, A, B, Q, P0{r,i});
        
        % extract upper triangular entries in a vector
        P0_iterates{r+1,i} = P0{r+1,i}(idx);
        err = norm(P0_iterates{r+1,i} - P0_iterates{r,i}, inf);
        r = r+1;
    end

    
end  

%% Cell resizing
P0 = P0(1:r);
P0_iterates = P0_iterates(1:r);

%% Checking contraction
for k = 2:r-1
    L_k = norm(P0{k+1} - P0{k}, 'fro') / norm(P0{k} - P0{k-1}, 'fro');
    if L_k >= 1
        disp('Warning: Not a contraction (L_k >= 1)');
    end
end

%% Plotting evolution of indices of the matrix
figure;
rows = ceil(sqrt(nsym));
cols = ceil(nsym / rows);

% Initialize a matrix to store all data (r rows x nsym columns)
all_data = zeros(r, nsym);

for i = 1:nsym
    subplot(rows, cols, i);

    for j = 1:m
        x = cellfun(@(v) v(i), P0_iterates(:, j)); % Extract values for plotting
        plot(1:r, x, '-o', 'LineWidth', 2);
        hold on;
    end

    % Store data for this i in the corresponding column
    all_data(:, i) = x;

    xlabel('iteration index $r$', 'FontSize', 10);
    ylabel(['$P_0(' num2str(i) ')$'], 'Interpreter', 'latex', 'FontSize', 10);
end

% Save all data to a single text file
filename = 'P0_iterates.txt';
fileID = fopen(filename, 'w');

% Write header
fprintf(fileID, 'Iteration ');
for i = 1:nsym
    fprintf(fileID, 'P0_%d ', i);
end
fprintf(fileID, '\n');

% Write data
for row = 1:r
    fprintf(fileID, '%d ', row); % Iteration index
    fprintf(fileID, '%.6f ', all_data(row, :)); % All P0 values for this iteration
    fprintf(fileID, '\n');
end

fclose(fileID);
%% Plotting Frobenius norms of all combinations (mC2)
% %Preallocate to store Frobenius norms
% frobenius_norms_P = zeros(max_iter, nchoosek(m, 2)); % Store norms for all pairs of P0 and P0 tilde 
% num_iter = max_iter;
% 
% % Iterate through iterations and calculate Frobenius norm for all pairs
% pair_idx = 1;
% for i = 1:m-1
%     for j = i+1:m
%         for k = 1:num_iter
%             % Calculate Frobenius norm between matrices from initial guesses i and j
%             frobenius_norms_P(k, pair_idx) = sqrt(trace((P0{k,i} - P0{k,j})'*(P0{k,i} - P0{k,j})));
%         end
%         pair_idx = pair_idx + 1;
%     end
% end
% 
% % Plot Frobenius norm trajectories
% figure;
% hold on;
% for i = 1:size(frobenius_norms_P, 2)
%     plot(1:num_iter, frobenius_norms_P(:, i));
% end
% hold off;
% title('Frobenius Norm Trajectories for $Y_0$');
% xlabel('Iteration $k$');
% ylabel('Frobenius Norm');

%% Forward propagating the Matrix ODEs
t_initial = 0; 
t_final = 1;

% integrator for Matrix Ricatti Equation P
N = 500; % number of grid points
t_fine = linspace(t_initial, t_final, N);

P_0 = zeros(nx);
for i = 1:m
    P_0 = P_0 + P0{end,i}; % initial condition for P
end
P_0 = P_0/m;

[t,P] = ode45(@(t,P) LyapODE_P(t,P,A,B,Q), t_fine, P_0(:));

% integrator for Matrix Ricatti Equation H
H0 = inv(Param.Sigma_0) - P_0; 
[t,H] = ode45(@(t,H) LyapODE_H(t,H,A,B,Q), t_fine, H0(:));

% Calculating evolution of Sigma matrices at time t
P_t = cell(length(t),1);
H_t = cell(length(t),1);
Sigma_t = cell(length(t),1);

for i = 1:length(t) 
    P_t{i} = reshape(P(i,:), size(A));
    H_t{i} = reshape(H(i,:), size(A));
    Sigma_t{i} = (H_t{i} + P_t{i})\eye(nx); % sigma at time t
end

%% Sampling initial points (DI)
k = 5; % number of initial points
x0 = -0.5 + 1*rand(k, 2);
dt = 0.0005;

u_t = zeros(1,N); % initialising control input
x_t = cell(k,1);

for i=1:k
   x_t{i} = Euler_Maruyama_SDE(x0(i,:)', dt, N, A, B, P_t); % Numerical approximation of x_t
end

%% Plotting matrices (2D-Ellipsoid)
figure(3)
q = [0;0]; % center of ellipsoid
RGB_row_vector_d = [0 1 0]; % for desired covariance matrix
RGB_row_vector_t = [0 1 1]; % for covariance matrix at time t

for i = 1:length(t) 
    S = Sigma_t{i}; 
    disp(sqrt(trace((S-Param.Sigma_d)'*(S-Param.Sigma_d)))/2);
    PlotEllipsoid2D(q, S, t(i), RGB_row_vector_t, 0.01);
    hold on;
    set(gca,'FontSize',30)
end
PlotEllipsoid2D(q, Param.Sigma_d, 1.001, RGB_row_vector_d, 0.2);
hold on;

% Noisy trajectory
for i=1:k
    plot3(t, x_t{i}(1,:), x_t{i}(2,:), 'LineWidth', 2); % plotting the noisy trajectory of x_t
end

grid on;
view(3);  % Set view to 3D perspective
xlabel('time','FontSize',30)
ylabel('Position x','FontSize',30)
zlabel('Position y','FontSize',30)

%% Sampling Points and Noisy Trajectories (CW)
% k = 5; % number of initial points
% %x0 = -0.5 + 1*rand(k, 6);
% x0 = -0.25 + 0.5*rand(k, 6);
% dt = 0.0005;
% 
% u_t = zeros(1,N); % initialising control input
% x_t = cell(k,1);
% 
% for i=1:k
%    x_t{i} = Euler_Maruyama_SDE(x0(i,:)', dt, N, A, B, P_t); % Numerical approximation of x_t
% end
% %% Plotting 3-D Ellipsoids (Position) (Example 4 - CW Equations)
% % Define a custom colormap from blue to yellow
% num_shades = 256; % Number of shades
% cmap = [linspace(0, 1, num_shades)', linspace(0, 1, num_shades)', linspace(1, 0, num_shades)']; 
% 
% figure(3)
% colormap(cmap); % Apply the custom colormap
% clim([0, 1]); % Normalize color scale (since t is between 0 and 1)
% colorbar; % Add colorbar to the figure
% hold on;
% 
% % For Position
% num_points = 10;  % Number of ellipsoids
% indices = sort([1, length(t), randi([2, length(t)-1], 1, num_points-2)]); 
% q = zeros(3, N); 
% j = 0; % Counter for translating the center
% translation_values = zeros(1, num_points); % Store y-axis translations for interpolation
% 
% % For storing Mean and Covariances of the ellipsoids
% q_pos = cell(num_points,1);
% S_pos = cell(num_points,1);
% 
% for idx = 1:num_points
%     i = indices(idx);
%     S = Sigma_t{i}; % Extract the covariance matrix Sigma(t)
%     S_pos{idx} = S(1:nx/2, 1:nx/2); % Extract the subblock corresponding to position
%     translation_values(idx) = 1.5 * j; % Compute translation along y-axis
%     q(:, i) = [0; translation_values(idx); 0]; % Center of ellipsoid
%     q_pos{idx} = q(:,i);
%     j = j + 1;
% 
%     % Map t_val to the colormap
%     t_val = t(i);
%     num_colors = size(cmap, 1);
%     color_idx = round(t_val * (num_colors - 1)) + 1; % Scale t_val to colormap indices
%     color = cmap(color_idx, :); % Get corresponding color from colormap
% 
%     PlotEllipsoid3D(q_pos{idx}, S_pos{idx}, color)
%     hold on;
%     set(gca, 'FontSize', 30)
% end
% 
% % ** Interpolate translation for all N points **
% translations = interp1(indices, translation_values, 1:N, 'linear', 'extrap');
% 
% % Apply translation to noisy trajectory
% qn = zeros(3, N);
% x_t_pshifted = cell(k,1);
% 
% for i = 1:N
%     qn(:, i) = [0; translations(i); 0]; % Apply interpolated translation
% 
%     for l = 1:k  % Loop over each cell
%         x_t_pshifted{l}(1:nx/2, i) = x_t{l}(1:nx/2, i) + qn(:, i);  % Apply translation
%     end
% end
% 
% % Plot final reference ellipsoid in red
% PlotEllipsoid3D(q(:, end), Param.Sigma_d(1:nx/2, 1:nx/2), [1, 0, 0]);
% grid on;
% hold on;
% 
% % Plot noisy trajectory
% for i = 1:k
%     plot3(x_t_pshifted{i}(1, :), x_t_pshifted{i}(2, :), x_t_pshifted{i}(3, :), 'LineWidth', 2); % Plot noisy trajectory
%     hold on;
% end
% 
% view(3);  % Set view to 3D perspective
% xlabel('Position x', 'FontSize', 30)
% ylabel('Position y', 'FontSize', 30)
% zlabel('Position z', 'FontSize', 30)
% 
% %% Plotting 3-D Ellipsoids (Velocity)
% % We'll only plot a certain number of ellipsoids.
% % Define a custom colormap from blue to yellow
% num_shades = 256; % Number of shades
% cmap = [linspace(0, 1, num_shades)', linspace(0, 1, num_shades)', linspace(1, 0, num_shades)']; 
% % This line creates an Nx3 matrix transitioning from blue [0 0 1] to yellow [1 1 0]
% 
% figure(4)
% colormap(cmap); % Apply the custom colormap
% clim([0, 1]); % Normalize color scale (since t is between 0 and 1)
% colorbar; % Add colorbar to the figure
% hold on;
% 
% % For storing Mean and Covariances of the ellipsoids
% S_vel = cell(num_points,1);
% j = 1; 
% 
% for i = indices  
%     S = Sigma_t{i}; % Extract the covariance matrix Sigma(t)
%     S_vel{j} = S(nx/2+1:nx, nx/2+1:nx); % Extract the subblock corresponding to position
% 
%     t_val = t(i); % time at which the ellipsoid is to be generated
% 
%     % Map t_val to the colormap
%     num_colors = size(cmap, 1);
%     color_idx = round(t_val * (num_colors - 1)) + 1; % Scale t_val to colormap indices
%     color = cmap(color_idx, :); % Get corresponding color from colormap
%     PlotEllipsoid3D(q_pos{j}, S_vel{j}, color)
%     hold on;
%     set(gca,'FontSize',30)
%     j = j+1;
% end
% 
% % Apply translation to noisy trajectory velocity
% x_t_vshifted = cell(k,1);
% 
% for i = 1:N
%     for l = 1:k  % Loop over each cell
%         x_t_vshifted{l}(1:nx/2, i) = x_t{l}(nx/2+1:end, i) + qn(:, i);  % Apply translation
%     end
% end
% 
% PlotEllipsoid3D(q(:,end), Param.Sigma_d(nx/2+1:nx, nx/2+1:nx), [1, 0, 0]);
% hold on;
% 
% % Plot noisy trajectory
% for i = 1:k
%     plot3(x_t_vshifted{i}(1, :), x_t_vshifted{i}(2, :), x_t_vshifted{i}(3, :), 'LineWidth', 2); % Plot noisy trajectory
%     hold on;
% end
% 
% grid on;
% view(3);  % Set view to 3D perspective
% xlabel('Velocity x','FontSize',30)
% ylabel('Velocity y','FontSize',30)
% zlabel('Velocity z','FontSize',30)

%% Optimal Control Plots (DI)
u = cell(k,1);
figure;

for v = 1:k
    for i = 1:length(t) 
        u{v}(:,i) = -B' * P_t{i} * x_t{v}(:,i);
    end

    % Plot u_x over time
    plot(t, u{v}(1,:), 'LineWidth', 2);
    hold on;
    xlabel('Time t');
    ylabel('u_x');
    title('Control Input u_x vs Time');
    grid on;

    % Save data to text file
    filename = sprintf('control_input_u_%d.txt', v); % Unique file per v
    fileID = fopen(filename, 'w');
    fprintf(fileID, 'Time u\n'); % Header
    fprintf(fileID, '%.6f %.6f\n', [t'; u{v}(1,:)]);
    fclose(fileID);
end

%% Optimal Control Plots (CW)
% u = cell(k,1);
% figure;
% 
% for v = 1:k
%     for i = 1:length(t) 
%         u{v}(:,i) = -B' * P_t{i} * x_t{v}(:,i);
%     end
% 
%     % Save data to text file
%     filename = sprintf('control_input_u_%d.txt', v); % Create unique filenames
%     fileID = fopen(filename, 'w');
%     fprintf(fileID, 'Time u_x u_y u_z\n'); % Header
%     fprintf(fileID, '%.6f %.6f %.6f %.6f\n', [t'; u{v}]); % Write data
%     fclose(fileID);
% 
%     % Subplot for u_1
%     subplot(3,1,1);
%     plot(t, u{v}(1,:), 'LineWidth', 2);
%     hold on;
%     xlabel('Time t');
%     ylabel('u_x');
%     title('Control Input u_x vs Time');
%     grid on;
% 
%     % Subplot for u_2
%     subplot(3,1,2);
%     plot(t, u{v}(2,:), 'LineWidth', 2);
%     hold on;
%     xlabel('Time t');
%     ylabel('u_y');
%     title('Control Input u_y vs Time');
%     grid on;
% 
%     % Subplot for u_3
%     subplot(3,1,3);
%     plot(t, u{v}(3,:), 'LineWidth', 2);
%     hold on;
%     xlabel('Time t');
%     ylabel('u_z');
%     title('Control Input u_z vs Time');
%     grid on;
% end

%% Evaluating Terminal Cost
Sigma_1 = Sigma_t{end};
phi_t = sqrt(trace((Sigma_1-Param.Sigma_d)*(Sigma_1-Param.Sigma_d)'))/2;

%% DI
% Store Noisy trajectory
for i = 1:k
    % Save data to text file
    filename = sprintf('noisy_trajectory_%d.txt', i); % Create unique filenames
    fileID = fopen(filename, 'w');
    fprintf(fileID, 'Time p_x p_y \n'); % Header
    fprintf(fileID, '%.6f %.6f %.6f \n', [t'; x_t{i}]); % Write data
    fclose(fileID);
end

%% Store Ellipse Data (DI)
Sig_d = Param.Sigma_d;
save('Covariance.mat', 'Sigma_t');
save('Time.mat', 't');
save('Desired_Covariance.mat', 'Sig_d')

% %% CW 
% % Store Noisy Position (Ellipsoid)
% 
% for i = 1:k
% 
%     % Save data to text file
%     filename = sprintf('noisy_position_%d.txt', i); % Create unique filenames
%     fileID = fopen(filename, 'w');
%     fprintf(fileID, 'Time p_x p_y p_z\n'); % Header
%     fprintf(fileID, '%.6f %.6f %.6f %.6f\n', [t'; x_t_pshifted{i}]); % Write data
%     fclose(fileID);
% end
% 
% % Store Noisy Velocity (Ellipsoid)
% for i = 1:k
% 
%     % Save data to text file
%     filename = sprintf('noisy_velocity_%d.txt', i); % Create unique filenames
%     fileID = fopen(filename, 'w');
%     fprintf(fileID, 'Time v_x v_y v_z\n'); % Header
%     fprintf(fileID, '%.6f %.6f %.6f %.6f\n', [t'; x_t_vshifted{i}]); % Write data
%     fclose(fileID);
% end
% 
% % Store Ellipsoid Data (CW)
% 
% save('Mean.mat', 'q_pos');
% save('Position_Covariance.mat', 'S_pos');
% save('Velocity_Covariance.mat', 'S_vel');