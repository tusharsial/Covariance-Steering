# Optimal Covariance Steering Algorithm in Continuous Time

## Overview

This repository contains the implementation of an **Optimal Covariance Steering Algorithm in Continuous Time with Hilbert-Schmidt Terminal Cost** for Linear Stochastic Systems over a finite time horizon. This work was developed under the supervision of **Professor Abhishek Halder**.

### Motivation

While there has been growing literature on fixed-horizon LQ covariance steering problems with terminal cost for discrete-time systems, the continuous-time version remains relatively unexplored. A key reason for this imbalance lies in computational tractability: the discrete-time formulation naturally leads to a semidefinite program solvable using off-the-shelf interior-point solvers, whereas the continuous-time formulation with terminal cost gives rise to a coupled nonlinear system of matricial ODEs, for which a principled and computationally efficient algorithm has remained unclear.

This research introduces a soft constraint via the **Hilbert-Schmidt Terminal Cost**, together with a quadratic cost function for control input and state. The necessary conditions of optimality lead to a coupled matrix ODE two-point boundary value problem with nonlinear split boundary conditions. To solve this system, we have designed a **Matricial Recursive Algorithm** with fast convergence rate, grounded in linear fractional transformations parameterized by the state transition matrix of the associated Hamiltonian system.

### Algorithm Validation

The proposed algorithm was tested and validated on a close-proximity rendezvous scenario by modeling the relative motion of a service spacecraft to a target satellite in LEO using **Clohessy-Wiltshire dynamics** with stochastic disturbances.

For complete details about the algorithm development, proof of convergence, and theoretical foundations, please refer to our manuscript available on [arXiv](YOUR_ARXIV_LINK_HERE).

---

## Requirements

- **MATLAB** (no additional toolboxes required)
- **Python** (optional, for generating publication-ready plots)
- **Hardware**: No high-performance computing required
  - *Note: The code was tested on [YOUR_CPU_MODEL], [NUMBER] cores, [GPU_MODEL if applicable], [RAM_AMOUNT]. Execution time was under 10 seconds.*

---

## Installation & Setup

### Step 1: Clone the Repository

Clone this repository to your preferred location:

```bash
git clone https://github.com/YOUR_USERNAME/YOUR_REPO_NAME.git
cd YOUR_REPO_NAME
```

### Step 2: Navigate to the MATLAB Code Directory

```bash
cd Code_base/Matlab_Codes_Covariance/
```

This directory contains all the required scripts to solve the optimal covariance steering problem.

---

## Running the Code

### Step 3: Execute the Main Script

Run the main MATLAB script:

```matlab
FTLQ_Covariance_Frobenius.m
```

This script solves the Fixed-time Linear Quadratic Covariance Steering problem in continuous time with a terminal cost measured in Hilbert-Schmidt (Frobenius) norm error between the desired and controlled terminal covariances.

### Step 4: Configure Your Example

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

This work lays the algorithmic groundwork toward handling more general terminal costs. Ongoing and future extensions include:

- **Wasserstein Terminal Distance**: Development of a custom algorithm to solve the boundary value problem with Wasserstein terminal cost
- **Generalized Noise Channels**: Extension to cases where noise does not enter through the same channel as the control input

The repository will be updated as these extensions are completed.

---

## Issues & Contributions

If you encounter any bugs or difficulties while executing the code, please open an issue or contact me directly. Your feedback is valuable and helps improve the implementation.

---

## Citation

If you use this code in your research, please cite our work:

```bibtex
@article{YOUR_CITATION_KEY,
  title={YOUR_PAPER_TITLE},
  author={YOUR_NAME and Abhishek Halder},
  journal={arXiv preprint arXiv:XXXX.XXXXX},
  year={2024/2025}
}
```

---

## License

[Add your license here, e.g., MIT, GPL, etc.]

---

## Contact

[Your Name]  
[Your Email]  
[Your Website/LinkedIn]

**Advisor**: Professor Abhishek Halder
