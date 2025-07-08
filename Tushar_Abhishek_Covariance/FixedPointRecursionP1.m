function P1 = FixedPointRecursionP1(H1, Sig_d, P1, nx)
%P1 = -inv(2*Sig_d)*(P1*P1 + P1*H1 + 2*(Sig_d*H1 - eye(nx)));
P1 = inv(-P1/2 + Sig_d) - H1;