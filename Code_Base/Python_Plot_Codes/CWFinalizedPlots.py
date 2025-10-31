import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from os.path import dirname, join as pjoin
import scipy.io
import GenerateRandomSPDmatrix
import PlotEllipsoid

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


#====================================================
# P0 recursion convergence plot
#====================================================

#fig, ax = plt.subplots()

# IterNum, P01, P02, P03 = np.loadtxt('P0_iterates.txt', unpack=True, dtype=float, usecols=range(4),skiprows=1)

# for k in range(21):
#     IterNum, Pcomponent = np.loadtxt('P0_iterates.txt', unpack=True, dtype=float, usecols=(0,k+1),skiprows=1)
#     plt.plot(IterNum, Pcomponent, 'k', lw=1, markerfacecolor='none',alpha=0.5)
#     plt.plot(IterNum[0], Pcomponent[0], 'ko', lw=1, markerfacecolor='none',alpha=0.5)

# plt.plot(IterNum, P01, 'r-o', lw=2, markerfacecolor='none')
# plt.plot(IterNum, P02, 'b-s', lw=2, markerfacecolor='none')
# plt.plot(IterNum, P03, 'g-d', lw=2, markerfacecolor='none')

# ax.set_ylabel(r"Iterates of $\bm{P}_{0}$")
# ax.set_xlabel(r"Recursion index")

# ax.autoscale (enable=True, axis='both', tight=1)

#ax.legend([r"$\bm{P}_{0}(1,1)$", r"$\bm{P}_{0}(1,2)$", r"$\bm{P}_{0}(2,2)$"], frameon=False)

# ax.spines['top'].set_visible(False)
# ax.spines['right'].set_visible(False)

# plt.savefig('Recursion_CW.png', dpi=300)

# =================================================================
# # # INPUT SAMPLE PATHS
# =================================================================
# t, u11, u12, u13 = np.loadtxt('control_input_u_1.txt', unpack=True, dtype=float, usecols=(0,1,2,3),skiprows=1)
# u21, u22, u23 = np.loadtxt('control_input_u_2.txt', unpack=True, dtype=float, usecols=(1,2,3),skiprows=1)
# u31, u32, u33 = np.loadtxt('control_input_u_3.txt', unpack=True, dtype=float, usecols=(1,2,3),skiprows=1)
# u41, u42, u43 = np.loadtxt('control_input_u_4.txt', unpack=True, dtype=float, usecols=(1,2,3),skiprows=1)
# u51, u52, u53 = np.loadtxt('control_input_u_5.txt', unpack=True, dtype=float, usecols=(1,2,3),skiprows=1)

# ax1 = plt.subplot(311)
# plt.plot(t,u11, color='#1f77b4', lw=3)
# plt.plot(t, u21, color='#2ca02c', lw=3)
# plt.plot(t, u31, color='#9467bd', lw=3)
# plt.plot(t, u41, color='#c7c7c7', lw=3)
# plt.plot(t, u51, color='#17becf', lw=3)
# plt.tick_params('x',labelbottom=False)
# ax1.set_ylabel(r"$u_{1t}^{\mathrm{opt}}$",rotation="horizontal")
# ax1.set_yticks([-0.5, 0, 0.5])

# ax2 = plt.subplot(312, sharex=ax1)
# plt.plot(t,u12, color='#1f77b4', lw=3)
# plt.plot(t, u22, color='#2ca02c', lw=3)
# plt.plot(t, u32, color='#9467bd', lw=3)
# plt.plot(t, u42, color='#c7c7c7', lw=3)
# plt.plot(t, u52, color='#17becf', lw=3)
# plt.tick_params('x',labelbottom=False)
# ax2.set_ylabel(r"$u_{2t}^{\mathrm{opt}}$",rotation="horizontal")
# ax2.set_yticks([-0.5, 0, 0.5])

# ax3 = plt.subplot(313, sharex=ax1)
# plt.plot(t,u13, color='#1f77b4', lw=3)
# plt.plot(t, u23, color='#2ca02c', lw=3)
# plt.plot(t, u33, color='#9467bd', lw=3)
# plt.plot(t, u43, color='#c7c7c7', lw=3)
# plt.plot(t, u53, color='#17becf', lw=3)
# plt.tick_params('x',labelbottom=False)
# ax3.set_ylabel(r"$u_{3t}^{\mathrm{opt}}$",rotation="horizontal")
# ax3.set_yticks([-0.5, 0, 0.5])

# ax3.set_xlabel(r"$t$")
# ax3.set_xticks([0, 1])
# ax3.set_xticklabels([0, 1])

# ax1.spines['top'].set_visible(False)
# ax1.spines['right'].set_visible(False)
# ax1.spines['bottom'].set_visible(False)
# ax2.spines['top'].set_visible(False)
# ax2.spines['right'].set_visible(False)
# ax2.spines['bottom'].set_visible(False)
# ax3.spines['top'].set_visible(False)
# ax3.spines['right'].set_visible(False)
# ax3.spines['bottom'].set_visible(False)

# # ax1.autoscale (enable=True, axis='both', tight=1)
# # ax2.autoscale (enable=True, axis='both', tight=1)
# # ax3.autoscale (enable=True, axis='both', tight=1)


# plt.savefig('Control_Effort_CW.png', dpi=300)

# # ==============================================================
# # STATE SAMPLE PATHS AND COVARIANCES
# # ==============================================================
# x11, x12, x13 = np.loadtxt('noisy_position_1.txt', unpack=True, dtype=float, usecols=(1,2,3),skiprows=1)
# x21, x22, x23 = np.loadtxt('noisy_position_2.txt', unpack=True, dtype=float, usecols=(1,2,3),skiprows=1)
# x31, x32, x33 = np.loadtxt('noisy_position_3.txt', unpack=True, dtype=float, usecols=(1,2,3),skiprows=1)
# x41, x42, x43 = np.loadtxt('noisy_position_4.txt', unpack=True, dtype=float, usecols=(1,2,3),skiprows=1)
# x51, x52, x53 = np.loadtxt('noisy_position_5.txt', unpack=True, dtype=float, usecols=(1,2,3),skiprows=1)

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')

# n_Ellipsoids = 10
# dim = 3

# class structtype():
#     pass
# SetOfEllipsoids = [ structtype() for i in range(n_Ellipsoids)]

# CenterVector = scipy.io.loadmat('Mean.mat')
# CenterVector = CenterVector['q_pos']

# ShapeMatrix = scipy.io.loadmat('Position_Covariance.mat')
# ShapeMatrix = ShapeMatrix['S_pos']

# # plot controlled covariances
# for i in range(n_Ellipsoids):
#     #SetOfEllipsoids[i].Center  = np.ravel(np.random.rand(dim,1))
#     SetOfEllipsoids[i].Center  = np.vstack(CenterVector[i]).flatten()

#     #SetOfEllipsoids[i].Shape  = GenerateRandomSPDmatrix.GenRandSPD(dim)
#     SetOfEllipsoids[i].Shape  = np.vstack(ShapeMatrix[i])
#     # for 2D
#     #PlotEllipsoid.TwoDim(SetOfEllipsoids[i].Center,SetOfEllipsoids[i].Shape,ax)
#     # for 3D
#     PlotEllipsoid.ThreeDim(SetOfEllipsoids[i].Center,SetOfEllipsoids[i].Shape,ax,'k')
# # plot the desired terminal covariance
# SigmaDesired = [[1.6431, 1.1138, 1.5453],
#  [1.1138, 1.9581, 1.4418],
#  [1.5453, 1.4418, 3.9142]]

# PlotEllipsoid.ThreeDim(SetOfEllipsoids[-1].Center,SigmaDesired,ax,'r')

# ax.plot3D(x11, x12, x13, color='#1f77b4', lw=1.5)
# ax.plot3D(x11[0], x12[0], x13[0], 'o', markeredgecolor='#1f77b4', markerfacecolor='none', lw=1.5)
# ax.plot3D(x21, x22, x23, color='#2ca02c', lw=1.5)
# ax.plot3D(x21[0], x22[0], x23[0], 'o', markeredgecolor='#2ca02c', markerfacecolor='none',lw=1.5)
# ax.plot3D(x31, x32, x33, color='#9467bd', lw=1.5)
# ax.plot3D(x31[0], x32[0], x33[0], 'o', markeredgecolor='#9467bd', markerfacecolor='none',lw=1.5)
# ax.plot3D(x41, x42, x43, color='k', lw=1.5)
# ax.plot3D(x41[0], x42[0], x43[0], 'o', markeredgecolor='k', markerfacecolor='none',lw=1.5)
# ax.plot3D(x51, x52, x53, color='#17becf', lw=1.5)
# ax.plot3D(x51[0], x52[0], x53[0], 'o', markeredgecolor='#17becf', markerfacecolor='none',lw=1.5)

# ax.set_xlabel(r"$x_{1t}^{\mathrm{opt}}$ [m]")
# ax.set_ylabel(r"$x_{3t}^{\mathrm{opt}}$ [m]")
# ax.zaxis.set_rotate_label(False)
# ax.set_zlabel(r"$x_{5t}^{\mathrm{opt}}$ [m]",rotation=90)
# ax.set_xticks([-1, 0, 1])
# ax.set_yticks([0, 10, 15])
# ax.set_zticks([-0.5, 0, 0.5])
# ax.grid(False)

# ax.azim = 19
# ax.elev = 9

# # For debugging
# #print(SetOfEllipsoids[0].Shape)

# # add rotated texts
# ax.text(0.6,0.7,0.8,r"$t=0.1$", rotation=85, fontsize=12, rotation_mode='anchor',zdir=(1,1,0))
# ax.text(0.6,3.4,0.82,r"$t=0.2$", rotation=85, fontsize=12, rotation_mode='anchor',zdir=(1,1,0))
# ax.text(0.6,5.15,0.78,r"$t=0.3$", rotation=85, fontsize=12, rotation_mode='anchor',zdir=(1,1,0))
# ax.text(0.6,6.52,0.73,r"$t=0.4$", rotation=85, fontsize=12, rotation_mode='anchor',zdir=(1,1,0))
# ax.text(0.6,7.7,0.69,r"$t=0.5$", rotation=85, fontsize=12, rotation_mode='anchor',zdir=(1,1,0))
# ax.text(0.6,9.3,0.7,r"$t=0.6$", rotation=85, fontsize=12, rotation_mode='anchor',zdir=(1,1,0))
# ax.text(0.6,10.6,0.68,r"$t=0.7$", rotation=85, fontsize=12, rotation_mode='anchor',zdir=(1,1,0))
# ax.text(0.6,12.1,0.68,r"$t=0.8$", rotation=85, fontsize=12, rotation_mode='anchor',zdir=(1,1,0))
# ax.text(0.6,13.2,0.65,r"$t=0.9$", rotation=85, fontsize=12, rotation_mode='anchor',zdir=(1,1,0))
# ax.text(0.6,15,0.66,r"$t=1$", rotation=85, fontsize=12, rotation_mode='anchor',zdir=(1,1,0))

# plt.savefig('StateCovariancesAndSamplePaths_CW.png', dpi=300)

# ==============================================================
# VELOCITY SAMPLE PATHS AND COVARIANCES
# ==============================================================
x11, x12, x13 = np.loadtxt('noisy_velocity_1.txt', unpack=True, dtype=float, usecols=(1,2,3),skiprows=1)
x21, x22, x23 = np.loadtxt('noisy_velocity_2.txt', unpack=True, dtype=float, usecols=(1,2,3),skiprows=1)
x31, x32, x33 = np.loadtxt('noisy_velocity_3.txt', unpack=True, dtype=float, usecols=(1,2,3),skiprows=1)
x41, x42, x43 = np.loadtxt('noisy_velocity_4.txt', unpack=True, dtype=float, usecols=(1,2,3),skiprows=1)
x51, x52, x53 = np.loadtxt('noisy_velocity_5.txt', unpack=True, dtype=float, usecols=(1,2,3),skiprows=1)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

n_Ellipsoids = 10
dim = 3

class structtype():
    pass
SetOfEllipsoids = [ structtype() for i in range(n_Ellipsoids)]

CenterVector = scipy.io.loadmat('Mean.mat')
CenterVector = CenterVector['q_pos']

ShapeMatrix = scipy.io.loadmat('Velocity_Covariance.mat')
ShapeMatrix = ShapeMatrix['S_vel']

# plot controlled covariances
for i in range(n_Ellipsoids):
    #SetOfEllipsoids[i].Center  = np.ravel(np.random.rand(dim,1))
    SetOfEllipsoids[i].Center  = np.vstack(CenterVector[i]).flatten()

    #SetOfEllipsoids[i].Shape  = GenerateRandomSPDmatrix.GenRandSPD(dim)
    SetOfEllipsoids[i].Shape  = np.vstack(ShapeMatrix[i])
    # for 2D
    #PlotEllipsoid.TwoDim(SetOfEllipsoids[i].Center,SetOfEllipsoids[i].Shape,ax)
    # for 3D
    PlotEllipsoid.ThreeDim(SetOfEllipsoids[i].Center,SetOfEllipsoids[i].Shape,ax,'k')
# plot the desired terminal covariance
SigmaDesired = [[2.1027, 1.3448, 0.9645],
 [1.3448, 1.7077, 0.7830],
 [0.9645, 0.7830, 1.5008]]

PlotEllipsoid.ThreeDim(SetOfEllipsoids[-1].Center,SigmaDesired,ax,'r')

ax.plot3D(x11, x12, x13, color='#1f77b4', lw=1.5)
ax.plot3D(x11[0], x12[0], x13[0], 'o', markeredgecolor='#1f77b4', markerfacecolor='none', lw=1.5)
ax.plot3D(x21, x22, x23, color='#2ca02c', lw=1.5)
ax.plot3D(x21[0], x22[0], x23[0], 'o', markeredgecolor='#2ca02c', markerfacecolor='none',lw=1.5)
ax.plot3D(x31, x32, x33, color='#9467bd', lw=1.5)
ax.plot3D(x31[0], x32[0], x33[0], 'o', markeredgecolor='#9467bd', markerfacecolor='none',lw=1.5)
ax.plot3D(x41, x42, x43, color='k', lw=1.5)
ax.plot3D(x41[0], x42[0], x43[0], 'o', markeredgecolor='k', markerfacecolor='none',lw=1.5)
ax.plot3D(x51, x52, x53, color='#17becf', lw=1.5)
ax.plot3D(x51[0], x52[0], x53[0], 'o', markeredgecolor='#17becf', markerfacecolor='none',lw=1.5)

ax.set_xlabel(r"$x_{2t}^{\mathrm{opt}}$ [m/s]")
ax.set_ylabel(r"$x_{4t}^{\mathrm{opt}}$ [m/s]")
ax.zaxis.set_rotate_label(False)
ax.set_zlabel(r"$x_{6t}^{\mathrm{opt}}$ [m/s]",rotation=90)
ax.set_xticks([-1, 0, 1])
ax.set_yticks([0, 15])
ax.set_zticks([-1, 0, 1])
ax.grid(False)

ax.azim = 12
ax.elev = 19

# For debugging
#print(SetOfEllipsoids[0].Shape)

# add rotated texts
ax.text(0.6,1.8,0.9,r"$t=0.1$", rotation=85, fontsize=12, rotation_mode='anchor',zdir=(1,1,0))
ax.text(0.6,3.4,0.9,r"$t=0.2$", rotation=85, fontsize=12, rotation_mode='anchor',zdir=(1,1,0))
ax.text(0.6,4.8,0.98,r"$t=0.3$", rotation=85, fontsize=12, rotation_mode='anchor',zdir=(1,1,0))
ax.text(0.6,6.2,1.15,r"$t=0.4$", rotation=85, fontsize=12, rotation_mode='anchor',zdir=(1,1,0))
ax.text(0.6,7.7,1.3,r"$t=0.5$", rotation=85, fontsize=12, rotation_mode='anchor',zdir=(1,1,0))
ax.text(0.6,9.2,1.32,r"$t=0.6$", rotation=85, fontsize=12, rotation_mode='anchor',zdir=(1,1,0))
ax.text(0.6,10.6,1.33,r"$t=0.7$", rotation=85, fontsize=12, rotation_mode='anchor',zdir=(1,1,0))
ax.text(0.6,12,1.35,r"$t=0.8$", rotation=85, fontsize=12, rotation_mode='anchor',zdir=(1,1,0))
ax.text(0.6,13.45,1.34,r"$t=0.9$", rotation=85, fontsize=12, rotation_mode='anchor',zdir=(1,1,0))
ax.text(0.6,15,1.31,r"$t=1$", rotation=85, fontsize=12, rotation_mode='anchor',zdir=(1,1,0))

plt.savefig('VelocityCovariancesAndSamplePaths_CW.png', dpi=300)