function Y0 = FixedPointRecursionY(Param, Y0)

% unpack the parameters
Phi = Param.Phi; Sig_0 = Param.Sigma_0; Sig_d = Param.Sigma_d;

% number of states
nx = size(Phi,1)/2;

% extract blocks of Phi
Phi_11 = Phi(1:nx,1:nx); Phi_12 = Phi(1:nx,nx+1:end);
Phi_21 = Phi(nx+1:end,1:nx); Phi_22 = Phi(nx+1:end,nx+1:end);

%% Recursion in Y0
L0 = (1/2)*inv(Sig_0) - Phi_12\Phi_11;
Z0 = L0 - Y0;
Sig_1 = PosDefProj(Phi_12*(Z0*Sig_0*Z0 - (1/4)*inv(Sig_0))*Phi_12');
disp(eig(Sig_1))
P1 = 2*(Sig_1 - Sig_d); % For Frobenius distance
%P1 = Sig_1\logm(Sig_1*inv(Sig_d)); % For Fisher-Rao distance
%P1 = sqrtm(Sig_d)*sqrtm(Sig_d^(-1/2)*inv(Sig_1)*Sig_d^(-1/2))*sqrtm(Sig_d) - eye(nx); % For Wasserstein Distance
Y0 = (P1*Phi_12 - Phi_22)\(Phi_21 - P1*Phi_11);