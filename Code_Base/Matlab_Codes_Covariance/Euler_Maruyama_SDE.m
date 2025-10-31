%%% This function implements the Euler–Maruyama (EM) method for the 
% approximate numerical solution of a stochastic differential equation (SDE).

function ys = Euler_Maruyama_SDE(x0, dt, N, A, B, P_t)
y_init  = x0;  % Initial y condition 

[nx, ~] = size(x0);
% Initialize the probability distribution for our
% random variable with mean 0 and stdev of sqrt(dt)
pd = makedist('Normal',0,sqrt(dt));

%ys = zeros(2,N);     % 1xN Matrix of zeros
ys = zeros(nx,N); 

% Generate a 2D Wiener process increment
ys(:,1) = y_init;

for i = 2:N
    y = ys(:,i-1);
    dw = random(pd,nx/2,1);
    P = P_t{i-1};
    ys(:,i) = y + (A-B*B'*P)*y*dt + B*dw;
end