function Pdot = LyapODE_P(t, P, A, B, Q)

P = reshape(P,size(A)); % reshape from vector to matrix

Pdot = -(A'*P + P*A - P*(B*B')*P +Q); % matrix ricatti differential equation

Pdot  = Pdot(:); % reshape from matrix to vector