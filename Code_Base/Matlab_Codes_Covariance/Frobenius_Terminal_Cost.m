function dist_fro = Frobenius_Terminal_Cost(Sigma_0, Sigma_d)
dist_fro = trace((Sigma_0 - Sigma_d)^2)/2; 