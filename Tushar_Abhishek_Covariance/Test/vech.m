function x = vech(X) % input must be a symmetric matrix X

% lower trinagularize
X = tril(X);

m = tril(true(size(X)));

x = X(m); % return column vector