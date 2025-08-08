import numpy as np
import scipy.linalg

def GenRandSPD(dim):
	# generate a random matrix
	M = np.random.rand(dim,dim) 
	# symmetrize
	M = 0.5*(M + M.transpose())
	# make it pos def
	M = M + dim*np.eye(dim)

	return M