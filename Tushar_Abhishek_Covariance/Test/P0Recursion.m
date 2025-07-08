close all; clear; clc;

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% parameters

dim = 2;

Sigma0 = randn(dim); Sigma0 = Sigma0*Sigma0'; Sigma0 = Sigma0 + (dim*eye(dim));
%X = rand(dim); Sigma0 = diag(rand(dim,1)) + X*X'; % Always positive definite

Sigma_d = randn(dim); Sigma_d = Sigma_d*(Sigma_d'); Sigma_d = Sigma_d + (dim*eye(dim));

A = rand(dim); B = rand(dim,1); Q = eye(dim); 

if rank(ctrb(A,B)) == dim
    M = [A -B*B';
        -Q -A'];
    
    Phi = expm(M);
    
    Phi11 = Phi(1:dim,1:dim); 
    Phi12 = Phi(1:dim,dim+1:end);
    Phi21 = Phi(dim+1:end,1:dim); 
    Phi22 = Phi(dim+1:end,dim+1:end); 
    
    % fixed point recursion
    maxiter = 999;
    
    numInitGuess = 5;

    P0 = cell(maxiter+1,numInitGuess); vechP0 = zeros(dim*(dim+1)/2,maxiter+1,numInitGuess);
    
    figure % phase portrait in \mathbb{S}^2 for n=2 states
    for j=1:numInitGuess
        % initial guess
        P0{1,j} = randn(dim); P0{1,j} = P0{1,j}*(P0{1,j})';
        vechP0(:,1,j) = vech(P0{1,j});

        for k=1:maxiter
    
            H0 = inv(Sigma0) - P0{k,j};
    
            H1 = -(Phi11' - H0*Phi12') \ (Phi21' - H0*Phi22');
    
            P1 = -0.5*(H1 + Sigma_d) + sqrtm((0.5*(H1 - Sigma_d))^2 + eye(dim));
    
            P0{k+1,j} = (P1*Phi12 - Phi22) \ (Phi21 - P1*Phi11);
    
            vechP0(:,k+1,j) = vech(P0{k+1,j});
        end
        plot3(vechP0(1,:,j),vechP0(2,:,j),vechP0(3,:,j),'-ko')
        hold on
        plot3(vechP0(1,1,j),vechP0(2,1,j),vechP0(3,1,j),'ro','MarkerSize',10,'MarkerFaceColor','r')
        hold on
        plot3(vechP0(1,end,j),vechP0(2,end,j),vechP0(3,end,j),'gs','MarkerSize',15,'MarkerFaceColor','g')
        hold on
    end
    hold off
    %axis tight
    grid on
    box on
    xlabel('$P_{0}(1,1)$','fontsize',30,'interpreter','latex');
    ylabel('$P_{0}(1,2)$','fontsize',30,'interpreter','latex');
    zlabel('$P_{0}(2,2)$','fontsize',30,'interpreter','latex');
    set(findall(gcf,'type','line'),'linewidth',1)
    set(gca,'FontSize',30)
    
    figure % components versus iter index, one subplot for each initial guess
    for j=1:numInitGuess
        subplot(numInitGuess,1,j)
        for i=1:(dim*(dim+1)/2)
            semilogx(1:maxiter+1, vechP0(i,:,j), '-ko')
            hold on
        end
        axis tight
        hold off
        %ylabel('Iterates of $P_{0}$','fontsize',30,'interpreter','latex');
        set(findall(gcf,'type','line'),'linewidth',1)
        set(gca,'FontSize',30)
    end
    xlabel('iteration index','fontsize',30,'interpreter','latex');
    set(gca,'FontSize',30)
else
    disp('NOT controllable');
end 