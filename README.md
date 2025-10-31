# 🧭 Fixed-Time LQ Covariance Steering with Hilbert–Schmidt Terminal Cost

This repository contains the MATLAB implementation of the **Fixed-Time Linear Quadratic (LQ) Covariance Steering** algorithm in **continuous time** with a **Hilbert–Schmidt (Frobenius) norm terminal cost**.  
The algorithm was developed as part of my **Master’s thesis under Professor Abhishek Halder**, and the details of the algorithm, convergence analysis, and theoretical foundations are presented in our accompanying manuscript — available [here on arXiv](<https://arxiv.org/abs/2510.21944>).

---

## Overview

This work addresses the **continuous-time LQ covariance steering problem** — a stochastic optimal control problem where the goal is to steer the covariance of a linear stochastic system over a fixed time horizon. Instead of a hard constraint, we've introduced a soft constraint here. 

While the discrete-time formulation is well-studied and can be solved using techniques including convex semidefinite programming, the continuous-time case has remained challenging due to the nonlinear coupling in the matrix ODE two-point boundary value problem.

To address this challenge, we propose a **Matricial Recursive Algorithm** with a **fast convergence rate**, leveraging **linear fractional transformations (LFTs)** parameterized by the **state transition matrix** of the associated **Hamiltonian system**.

The algorithm has been validated on various LIT systems. However, in the paper and code, we discuss the performance of the algorithm for the following 2 examples:
1. **Double Integrator (2D)** system, and  
2. **Noisy Clohessy–Wiltshire (CW)** equations (6D) for rendezvous of a chaser spacecraft with a target satellite under stochastic disturbances.

Our broader goal is to extend this framework for **Wasserstein terminal cost** and for cases where **noise enters through a channel different from the control input**.

---

## 🧩 Repository Structure

Code_base/
│
├── Matlab_Codes_Covariance/
│ ├── FTLQ_Covariance_Frobenius.m ← Main script (run this first)
│ ├── [other supporting MATLAB files]
│
├── Python_Plot_Codes/ ← Optional post-processing scripts
│ ├── DoubleIntegratorFinalizedPlots.py
│ ├── CWFinalizedPlots.py
│
├── Numerical_Example_Data/ ← Auto-generated directory (data files)
└── Plots/ ← Auto-generated directory (figures)


---

## 🚀 How to Run the Code

### Step 1: Clone the Repository
```bash
git clone https://github.com/<your-username>/<repo-name>.git



