function P1_pos = PosDefProj(P1)
[V,D] = eig(P1);
diag_elements = diag(D);
diag_elements(diag_elements<0) = 0;
D = diag(diag_elements);
P1_pos = V*D*V';