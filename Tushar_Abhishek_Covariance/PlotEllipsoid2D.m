function PlotEllipsoid2D(q, Q, t, RGB_row_vector,alp)

% Get some x/y values on the unit circle, for ellipse drawing purpose
angle = linspace(0,2*pi,500)'; xy = [sin(angle) cos(angle)];
% Transform using Cholesky decomposition
XY = xy*chol(Q);
% RGB_row_vector = [0 0 1] is blue, [1 0 0] is red 
%fill(q(1)+XY(:,1), q(2)+XY(:,2), RGB_row_vector) 
fill3(t * ones(size(XY, 1), 1), q(1) + XY(:,1), q(2) + XY(:,2), RGB_row_vector, 'FaceAlpha', alp, 'EdgeColor', [1 0 1], 'EdgeAlpha', 0.2);

view(3); % Set to 3D view