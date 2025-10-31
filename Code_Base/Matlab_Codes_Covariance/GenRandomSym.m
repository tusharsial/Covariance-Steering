function Y = GenRandomSym(a,b,n)

Y = a + (b-a)*rand(n);

Y  = 0.5*(Y + Y'); % symmetrize

%Y = (1/n)*Y; % make spectral rad = matrix 2 norm < 1

