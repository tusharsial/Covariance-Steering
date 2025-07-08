function Sig_1 = FixedPointRecursionSigma(Param, Sig_1)

% unpack the parameters
Phi = Param.Phi; Sig_0 = Param.Sigma_0; Sig_d = Param.Sigma_d;

% number of states
nx = size(Phi,1)/2;

% extract blocks of Phi
Phi_11 = Phi(1:nx,1:nx); Phi_12 = Phi(1:nx,nx+1:end);
Phi_21 = Phi(nx+1:end,1:nx); Phi_22 = Phi(nx+1:end,nx+1:end);

% % Recursion equations (1)
% K = sqrtm(Sig_0)*inv(Phi_12);
% M = sqrtm(inv(Sig_0))*sqrtm((eye(nx)/4 + K*Sig_1*K'))*sqrtm(inv(Sig_0));
% Y0 = inv(Sig_0)/2 - Phi_12\Phi_11 - M;
% P1 = (Phi_21 + Phi_22*Y0)*inv(Phi_11 + Phi_12*Y0);
% P1_pos = PosDefProj(P1);
% Sig_1 = P1_pos + Sig_d;
% 
% %disp(trace((P1_pos - P1)'*(P1_pos - P1)));
% %disp(eig(Y0+Phi_12\Phi_11));
% disp(eig(P1_pos));

% Recursion equations (2)
P1 = 2*(Sig_1 - Sig_d); % For Frobenius distance
%P1 = Sig_1\logm(Sig_1*inv(Sig_d)); % For Fisher-Rao distance
%P1 = sqrtm(Sig_d)*sqrtm(Sig_d^(-1/2)*inv(Sig_1)*Sig_d^(-1/2))*sqrtm(Sig_d) - eye(nx);
Y0 = (P1*Phi_12 - Phi_22)\(Phi_21 - P1*Phi_11);
Z0 = inv(Sig_0)/2 - Phi_12\Phi_11 - Y0;
L = (Z0*Sig_0*Z0 - inv(Sig_0)/2);
%L_pos = PosDefProj(L);
%Sig_1 = Phi_12*L_pos*Phi_12';
Sig_1 = PosDefProj(Phi_12*L*Phi_12');

disp(eig(Sig_1))
%disp(eig(inv(Sig_0)/2 - sqrtm(inv(Sig_0))*sqrtm((eye(nx)/4 + sqrtm(Sig_0)*inv(Phi_12)*Sig_1*(sqrtm(Sig_0)*inv(Phi_12))'))*sqrtm(inv(Sig_0))));