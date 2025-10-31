import numpy as np
import scipy.linalg as la
from matplotlib import rc
import matplotlib.pyplot as plt
# from scipy.integrate import odeint
import GenerateRandomSPDmatrix
import PlotEllipsoid
# import LownerJohn

dim = 2
n_Ellipsoids = 2

class structtype():
    pass
SetOfEllipsoids = [ structtype() for i in range(n_Ellipsoids)]

# for 2D plot
fig, axs = plt.subplots(1,1,sharex=True, sharey=True) 
# # for 3D plot
# fig = plt.figure()
# axs = fig.add_subplot(111, projection='3d')

for i in range(n_Ellipsoids):
	SetOfEllipsoids[i].Center  = np.ravel(np.random.rand(dim,1))
	SetOfEllipsoids[i].Shape  = GenerateRandomSPDmatrix.GenRandSPD(dim)
	# for 2D
	PlotEllipsoid.TwoDim(SetOfEllipsoids[i].Center,SetOfEllipsoids[i].Shape,axs)
	# for 3D
	# PlotEllipsoid.ThreeDim(SetOfEllipsoids[i].Center,SetOfEllipsoids[i].Shape,axs)

# print LownerJohn.MVOE_UnionOfEllipsoids(SetOfEllipsoids)

R = (la.inv(SetOfEllipsoids[0].Shape))*SetOfEllipsoids[1].Shape

print la.eig(R)

plt.savefig('2Dellipse.png', dpi=300) # for 2D
# plt.savefig('3Dellipsoid.png', dpi=300) # for 3D
