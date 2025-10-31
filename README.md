# Fixed Horizon Linear Quadratic Covariance Steering in Continuous Time with Hilbert-Schmidt Terminal Cost

## Overview

This repository contains the implementation of a **Fixed Horizon Linear Quadratic Covariance Steering in Continuous Time with Hilbert-Schmidt Terminal Cost**.  The algorithm was developed as part of my **Master’s thesis under Professor Abhishek Halder**, and the details of the algorithm, convergence analysis, and theoretical foundations are presented in our accompanying manuscript — available [here on arXiv](<https://arxiv.org/abs/2510.21944>).

---

### Motivation

While there has been a growing body of literature on fixed-horizon LQ covariance steering problems with terminal cost for discrete-time systems, the continuous-time version remains relatively unexplored. A key reason for this imbalance lies in computational tractability: the discrete-time formulation naturally leads to a semidefinite program that can be solved using off-the-shelf interior-point solvers, whereas the continuous-time formulation with terminal cost gives rise to a coupled nonlinear system of matrix ODEs, for which a principled and computationally efficient algorithm has remained unclear.

To address this challenge, we propose a **Matricial Recursive Algorithm** with a **fast convergence rate**, leveraging **linear fractional transformations (LFTs)** parameterized by the **state transition matrix** of the associated **Hamiltonian system**.

This research introduces a soft constraint via the **Hilbert-Schmidt Terminal Cost**, together with a quadratic cost function for control input and state. The necessary conditions of optimality lead to a coupled matrix ODE two-point boundary value problem with nonlinear split boundary conditions. To solve this system, we have designed a **Matricial Recursive Algorithm** with a fast convergence rate, grounded in linear fractional transformations parameterized by the state transition matrix of the associated Hamiltonian system.

### Algorithm Validation

The proposed algorithm was tested and validated on multiple systems. Below are 2 such examples:
1. **Double Integrator (2D)** system, and  
2. **Noisy Clohessy–Wiltshire (CW)** equations (6D) for spacecraft rendezvous under stochastic disturbances.

Our broader goal is to extend this framework to **Wasserstein terminal costs** and cases where **noise enters through a channel different from the control input**.

---

## Requirements

- **MATLAB** (no additional toolboxes required)
- **Python** (optional, for generating publication-ready plots)
- **Hardware**: No high-performance computing required
  - *Note: The code was tested on [AMD Ryzen 5-4600H (8 cores, 3.00 GHz)], [NUMBER] cores, [NVIDIA 1650 Ti 1600], [16 GB RAM]. Execution time was under 10 seconds.*

---

## Installation & Setup

### Step 1: Clone the Repository

Clone this repository to your preferred location:

### Step 2: Navigate to the MATLAB Code Directory

```bash
cd Code_base/Matlab_Codes_Covariance/
```

This directory contains all the required scripts to solve the optimal covariance steering problem.

---

## Running the Code

### Step 3: Configure Your Example

**IMPORTANT**: Before executing the code, please configure the following parameters in `FTLQ_Covariance_Frobenius.m`:

The script contains two numerical examples:
- **DI**: Double Integrator (2D)
- **Noisy CW**: Rendezvous scenario with Noisy Clohessy-Wiltshire equations (6D)

#### Configuration Parameters:

1. **Set Example (for plots)**
   - `Example = 'DI'` for Double Integrator
   - `Example = 'CW'` for Noisy CW

2. **Select Directory for data storage**
   - `current_dir = di_dir` for Double Integrator
   - `current_dir = noisy_cw_dir` for Noisy CW

3. **Set System Dimensions**
   - `nx = 2` for Double Integrator
   - `nx = 6` for Noisy CW

4. **Uncomment/Comment Sections**
   - **To run the "DI" example**: Uncomment all lines in the "DI" section and comment out all lines in the "Noisy CW" section
   - **To run the "Noisy CW" example**: Uncomment all lines in the "Noisy CW" section and comment out all lines in the "DI" section

### Step 4: Execute the Main Script

Run the main MATLAB script:

```matlab
FTLQ_Covariance_Frobenius.m
```

This script solves the Fixed-time Linear Quadratic Covariance Steering problem in continuous time with a terminal cost measured in Hilbert-Schmidt (Frobenius) norm error between the desired and controlled terminal covariances.

### Step 5: Output

The MATLAB script will automatically:
- Export numerical data as `.txt` files to `Numerical_Example_Data/` directory
- Save plots as `.png` files to `Plots/` directory
- Create these directories automatically if they don't exist

---

## Optional: Publication-Ready Plots (Python)

For generating publication-quality plots, Python scripts are provided in:

```bash
cd Code_base/Python_Plot_Codes/
```

**Scripts:**
- `DoubleIntegratorFinalizedPlots.py` - For DI example
- `CWFinalizedPlots.py` - For Noisy CW example

These scripts use the exported data from `Numerical_Example_Data/` as input to generate enhanced visualizations.

---

## Future Work

This work lays the algorithmic groundwork toward handling the Wasserstein terminal cost. Future extensions include:

- **Wasserstein Terminal Distance**: Development of a custom algorithm to solve the boundary value problem with Wasserstein terminal cost
- **Generalized Noise Channels**: Extension to cases where noise does not enter through the same channel as the control input

---

## Issues & Contributions

If you encounter any bugs or difficulties while executing the code, please open an issue or contact me directly. Your feedback is valuable and helps improve the implementation.

---

## Citation

If you use this code in your research, please cite our work:

```bibtex
@misc{sial2025fixedhorizonlinearquadratic,
      title={Fixed Horizon Linear Quadratic Covariance Steering in Continuous Time with Hilbert-Schmidt Terminal Cost}, 
      author={Tushar Sial and Abhishek Halder},
      year={2025},
      eprint={2510.21944},
      archivePrefix={arXiv},
      primaryClass={math.OC},
      url={https://arxiv.org/abs/2510.21944}, 
}
```

---

## Contact

[Tushar Sial]  
[tsial@iastate.edu]  
[https://tusharsial.github.io/]

**Advisor**: Professor Abhishek Halder
