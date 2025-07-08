function X = GenRandomPosDef(n)

X = rand(n); 
X = diag(rand(n,1)) + X*X';