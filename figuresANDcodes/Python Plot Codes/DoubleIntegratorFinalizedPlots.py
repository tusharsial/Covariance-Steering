import numpy as np
import scipy.linalg
import matplotlib.pyplot as plt
from os.path import dirname, join as pjoin
import scipy.io
import GenerateRandomSPDmatrix
import PlotEllipsoid

#====================================================
# Make plots beautiful
#====================================================

pts_per_inch = 72.27
# write "\the\textwidth" (or "\showthe\columnwidth" for a 2 collumn text)
text_width_in_pts = 300.0
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


#====================================================
# P0 recursion plot
#====================================================

#fig, ax = plt.subplots()

# IterNum, P01, P02, P03 = np.loadtxt('P0_iterates.txt', unpack=True, dtype=float, usecols=range(4))

# plt.plot(IterNum, P01, 'r-o', lw=2, markerfacecolor='none')
# plt.plot(IterNum, P02, 'b-s', lw=2, markerfacecolor='none')
# plt.plot(IterNum, P03, 'g-d', lw=2, markerfacecolor='none')

# ax.set_ylabel(r"Iterates of $\bm{P}_{0}$")
# ax.set_xlabel(r"Recursion index")

# # ax.autoscale (enable=True, axis='both', tight=1)

# ax.legend([r"$\bm{P}_{0}(1,1)$", r"$\bm{P}_{0}(1,2)$", r"$\bm{P}_{0}(2,2)$"], frameon=False)

# ax.spines['top'].set_visible(False)
# ax.spines['right'].set_visible(False)

# plt.savefig('DIplotPrecursion.png', dpi=300)


# # INPUT SAMPLE PATHS

# t, u1 = np.loadtxt('control_input_u_1.txt', unpack=True, dtype=float, usecols=range(2),skiprows=1)
# u2 = np.loadtxt('control_input_u_2.txt', unpack=True, dtype=float, usecols=1,skiprows=1)
# u3 = np.loadtxt('control_input_u_3.txt', unpack=True, dtype=float, usecols=1,skiprows=1)
# u4 = np.loadtxt('control_input_u_4.txt', unpack=True, dtype=float, usecols=1,skiprows=1)
# u5 = np.loadtxt('control_input_u_5.txt', unpack=True, dtype=float, usecols=1,skiprows=1)

# ax.set_prop_cycle(color=[
#     '#1f77b4', '#2ca02c', '#9467bd', '#c7c7c7', '#17becf'])
#     #bcbd22', '#dbdb8d',
#     # '#17becf', '#9edae5'])

# plt.plot(t, u1, color='#1f77b4', lw=2)
# plt.plot(t, u2, color='#2ca02c', lw=2)
# plt.plot(t, u3, color='#9467bd', lw=2)
# plt.plot(t, u4, color='#c7c7c7', lw=2)
# plt.plot(t, u5, color='#17becf', lw=2)

# ax.set_ylabel(r"$u_{t}^{\mathrm{opt}}$",rotation="horizontal")
# ax.set_xlabel(r"$t$")

# ax.spines['top'].set_visible(False)
# ax.spines['right'].set_visible(False)

# plt.savefig('DIplotInputPaths.png', dpi=300)


# ==============================================================
# STATE SAMPLE PATHS AND COVARIANCES
# ==============================================================
t, x11, x12 = np.loadtxt('noisy_trajectory_1.txt', unpack=True, dtype=float, usecols=(0,1,2),skiprows=1)
x21, x22 = np.loadtxt('noisy_trajectory_2.txt', unpack=True, dtype=float, usecols=(1,2),skiprows=1)
x31, x32 = np.loadtxt('noisy_trajectory_3.txt', unpack=True, dtype=float, usecols=(1,2),skiprows=1)
x41, x42 = np.loadtxt('noisy_trajectory_4.txt', unpack=True, dtype=float, usecols=(1,2),skiprows=1)
x51, x52 = np.loadtxt('noisy_trajectory_5.txt', unpack=True, dtype=float, usecols=(1,2),skiprows=1)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

n_Ellipsoids = 500
dim = 2

class structtype():
    pass
SetOfEllipsoids = [ structtype() for i in range(n_Ellipsoids)]

ShapeMatrix = scipy.io.loadmat('Covariance.mat')
ShapeMatrix = ShapeMatrix['Sigma_t']

angle = np.linspace(0, 2*np.pi, 500)

xy = np.array([[np.sin(angle)], [np.cos(angle)]]).transpose()

# plot controlled covariances
for i in range(n_Ellipsoids):
    SetOfEllipsoids[i].Center  = np.ravel(np.zeros(dim))
    # SetOfEllipsoids[i].Center  = np.vstack(CenterVector[i]).flatten()

    # SetOfEllipsoids[i].Shape  = GenerateRandomSPDmatrix.GenRandSPD(dim)
    SetOfEllipsoids[i].Shape  = np.vstack(ShapeMatrix[i])

    XY = xy.dot(np.linalg.cholesky(SetOfEllipsoids[i].Shape).transpose())
    # for 2D
    ax.plot3D(t[i]*np.ones(500), XY[:,0,0], XY[:,0,1], '-k', lw=1, alpha=0.1)
    # PlotEllipsoid.TwoDim(SetOfEllipsoids[i].Center,SetOfEllipsoids[i].Shape,ax,'none','k',0.5)
    # for 3D
    #PlotEllipsoid.ThreeDim(SetOfEllipsoids[i].Center,SetOfEllipsoids[i].Shape,ax,'k')
# plot the desired terminal covariance
SigmaDesired = scipy.io.loadmat('Desired_Covariance.mat')
SigmaDesired = SigmaDesired['Sig_d']
XYDesired = xy.dot(np.linalg.cholesky(SigmaDesired).transpose())
ax.plot3D(t[-1]*np.ones(500), XYDesired[:,0,0], XYDesired[:,0,1], '-r', lw=1)
# SigmaDesired = [[2.1027, 1.3448, 0.9645],
#  [1.3448, 1.7077, 0.7830],
#  [0.9645, 0.7830, 1.5008]]

# PlotEllipsoid.ThreeDim(SetOfEllipsoids[-1].Center,SigmaDesired,ax,'r')

ax.plot3D(t, x11, x12, color='#1f77b4', lw=1.5)
ax.plot3D(t[0], x11[0], x12[0], 'o', markeredgecolor='#1f77b4', markerfacecolor='none', lw=1.5)
ax.plot3D(t, x21, x22, color='#2ca02c', lw=1.5)
ax.plot3D(t[0], x21[0], x22[0], 'o', markeredgecolor='#2ca02c', markerfacecolor='none',lw=1.5)
ax.plot3D(t, x31, x32, color='#9467bd', lw=1.5)
ax.plot3D(t[0], x31[0], x32[0], 'o', markeredgecolor='#9467bd', markerfacecolor='none',lw=1.5)
ax.plot3D(t, x41, x42, color='k', lw=1.5)
ax.plot3D(t[0], x41[0], x42[0], 'o', markeredgecolor='k', markerfacecolor='none',lw=1.5)
ax.plot3D(t, x51, x52, color='#17becf', lw=1.5)
ax.plot3D(t[0], x51[0], x52[0], 'o', markeredgecolor='#17becf', markerfacecolor='none',lw=1.5)

ax.set_xlabel(r"$t$")
ax.set_ylabel(r"$x_{1t}^{\mathrm{opt}}$")
ax.zaxis.set_rotate_label(False)
ax.set_zlabel(r"$x_{2t}^{\mathrm{opt}}$",rotation=90)
ax.set_xticks([0, 1])
ax.set_yticks([-2, 2])
ax.set_zticks([-2, 2])
ax.grid(False)

ax.azim = 240
ax.elev = 19

# For debugging
# print(XY[:,0,0].size)

# add rotated texts
# ax.text(0.6,1.8,0.9,r"$t=0.1$", rotation=85, fontsize=12, rotation_mode='anchor',zdir=(1,1,0))
# ax.text(0.6,3.4,0.9,r"$t=0.2$", rotation=85, fontsize=12, rotation_mode='anchor',zdir=(1,1,0))
# ax.text(0.6,4.8,0.98,r"$t=0.3$", rotation=85, fontsize=12, rotation_mode='anchor',zdir=(1,1,0))
# ax.text(0.6,6.2,1.15,r"$t=0.4$", rotation=85, fontsize=12, rotation_mode='anchor',zdir=(1,1,0))
# ax.text(0.6,7.7,1.3,r"$t=0.5$", rotation=85, fontsize=12, rotation_mode='anchor',zdir=(1,1,0))
# ax.text(0.6,9.2,1.32,r"$t=0.6$", rotation=85, fontsize=12, rotation_mode='anchor',zdir=(1,1,0))
# ax.text(0.6,10.6,1.33,r"$t=0.7$", rotation=85, fontsize=12, rotation_mode='anchor',zdir=(1,1,0))
# ax.text(0.6,12,1.35,r"$t=0.8$", rotation=85, fontsize=12, rotation_mode='anchor',zdir=(1,1,0))
# ax.text(0.6,13.45,1.34,r"$t=0.9$", rotation=85, fontsize=12, rotation_mode='anchor',zdir=(1,1,0))
# ax.text(0.6,15,1.31,r"$t=1$", rotation=85, fontsize=12, rotation_mode='anchor',zdir=(1,1,0))


plt.savefig('StateCovariancesAndSamplePaths_DI.png', dpi=300)
