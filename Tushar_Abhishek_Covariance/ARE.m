function X = ARE(H_1, Sig_d)
nx = size(Sig_d);

I = eye(nx);
A = -(H_1 + Sig_d)/2;
B = I;
R = I;
E = I;
Q = I - (Sig_d*H_1 + H_1*Sig_d)/2;

[X, K, L] = icare(A, B, Q, R, [], E, []); % gives positive definite solution

%X = icare(-A, [], -Q, [], [], [], I); % gives negative definite solution
%(have to reject it but how

%% Test
%disp(rank(ctrb(A,B)));
%L = (H_1 - Sig_d)*X + X*(Sig_d - H_1) + H_1*Sig_d - Sig_d*H_1;
Sol = -(H_1+Sig_d)/2 + sqrtm(((H_1-Sig_d)/2)^2 + I);
%disp(L);
%fprintf('Eigenvalues of P1\n');
%disp(eig(X))