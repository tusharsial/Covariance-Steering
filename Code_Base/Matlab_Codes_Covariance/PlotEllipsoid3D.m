%%% Function for plotting 3D Covariance Ellipsoids

function PlotEllipsoid3D(q, Q, color)

N = 30;

% perform svd
[U,D,V] = svd(Q);

% get lengths of semi-axes
a = 1./diag(sqrt(D));

% generate centered ellipsoid
[X,Y,Z] = ellipsoid(0,0,0,a(1),a(2),a(3),N);

% rotate and translate
for k = 1:length(X)
    for j = 1:length(X)
            generic_pt = [X(k,j) Y(k,j) Z(k,j)]';
            new_pt = V * generic_pt;
            XX(k,j) = new_pt(1) + q(1);
            YY(k,j) = new_pt(2) + q(2);
            ZZ(k,j) = new_pt(3) + q(3);
    end
end

%% Plotting 
ellipsoid_plot = surf(XX, YY, ZZ, 'FaceColor', color, 'EdgeColor', 'none');

%shading interp
alpha(ellipsoid_plot,0.3);
axis tight
