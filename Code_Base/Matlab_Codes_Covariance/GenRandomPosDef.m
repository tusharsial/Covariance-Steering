% Function for generating an n x n random positive definite matrix 

function X = GenRandomPosDef(n)

X = rand(n); 
X = diag(rand(n,1)) + X*X';