import numpy as np
import scipy.linalg
from matplotlib import rc
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
#====================================================
# Make plots beautiful
#====================================================

pts_per_inch = 72.27
# write "\the\textwidth" (or "\showthe\columnwidth" for a 2 collumn text)
text_width_in_pts = 250.0
# inside a figure environment in latex, the result will be on the
# dvi/pdf next to the figure. See url above.
text_width_in_inches = text_width_in_pts / pts_per_inch
# figure.png or figure.eps will be intentionally larger, because it is prettier
inverse_latex_scale = 2
fig_proportion = (3.0 / 3.0)
csize = inverse_latex_scale * fig_proportion * text_width_in_inches
# always 1.0 on the first argument
fig_size = (1.0 * csize, 0.85 * csize)
# find out the fontsize of your latex text, and put it here
text_size = inverse_latex_scale * 9
label_size = inverse_latex_scale * 10
tick_size = inverse_latex_scale * 8
# learn how to configure:
# http://matplotlib.sourceforge.net/users/customizing.html
params = {'backend': 'ps',
          'axes.labelsize': 16,
          'legend.fontsize': tick_size,
          'legend.handlelength': 2.5,
          'legend.borderaxespad': 0,
          'axes.labelsize': label_size,
          'xtick.labelsize': tick_size,
          'ytick.labelsize': tick_size,
          'font.family': 'serif',
          'font.size': text_size,
          'font.serif': ['Computer Modern Roman'],
          'ps.usedistiller': 'xpdf',
          'text.usetex': True,
          'figure.figsize': fig_size,
          # # include here any neede package for latex
          # r"\usepackage{bm} \usepackage{amsmath}",
          }
plt.rcParams.update(params)
plt.rcParams['text.latex.preamble'] = r"\usepackage{bm} \usepackage{amsmath}"
fig = plt.figure(1, figsize=fig_size)  # figsize accepts only inches.
fig.subplots_adjust(left=0.04, right=0.98, top=0.93, bottom=0.15,
                    hspace=0.05, wspace=0.02)



def TwoDim(q,Q,axs,fcolor,ecolor,trnsprncy):

	angle = np.linspace(0, 2*np.pi, 500)

	xy = np.array([[np.sin(angle)], [np.cos(angle)]]).transpose()

	XY = xy.dot(np.linalg.cholesky(Q).transpose())

	# fig, axs = plt.subplots(1,1,sharex=True, sharey=True)
	# axs.fill(q[0]+XY[:,0,0], q[1]+XY[:,0,1], alpha=trnsprncy, facecolor=fcolor, edgecolor=ecolor, linewidth=1, zorder=1)

def ThreeDim(q,Q,axs,mycolor):

	# find the rotation matrix (V) and length of semi-axes (radii) from SVD
	U, D, V = np.linalg.svd(Q)
	radii = 1.0/np.sqrt(D)

	u = np.linspace(0.0, 2.0 * np.pi, 200)
	v = np.linspace(0.0, np.pi, 200)
	x = radii[0] * np.outer(np.cos(u), np.sin(v))
	y = radii[1] * np.outer(np.sin(u), np.sin(v))
	z = radii[2] * np.outer(np.ones_like(u), np.cos(v))
	
	for i in range(len(x)):
		for j in range(len(x)):
			[x[i,j],y[i,j],z[i,j]] = np.dot([x[i,j],y[i,j],z[i,j]], V) + q

	axs.plot_wireframe(x, y, z,  rstride=4, cstride=4, color=mycolor, alpha=0.05)	
