function Hdot = LyapODE_H(t, H, A, B, Q)

H = reshape(H,size(A)); % reshape from vector to matrix

Hdot = -(A'*H + H*A + H*(B*B')*H - Q); % matrix ricatti differential equation

Hdot  = Hdot(:); % reshape from matrix to vector