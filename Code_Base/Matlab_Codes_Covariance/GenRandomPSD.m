% Function for generating an n x n random positive semi-definite matrix.

function X = GenRandomPSD(n)
X = randn(n);

% Compute A = B'*B (which is always PSD)
X = X'*X;
end