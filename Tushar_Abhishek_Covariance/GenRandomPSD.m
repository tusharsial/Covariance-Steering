function X = GenRandomPSD(n)

X = randn(n);

% Compute A = B'*B (which is always PSD)
X = X'*X;
end