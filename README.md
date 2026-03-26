# Optimization Framework using AMPL

This repository provides the AMPL implementation of the optimization framework presented in the manuscript:
"Optimizing genetic gain and diversity in seed orchards using an algebraic modeling language for mathematical programming."

The framework formulates seed orchard deployment as a mixed-integer quadratic programming (MIQP) problem to maximize genetic gain while maintaining constraints on genetic diversity.

---

## Overview
- Algebraic and transparent model structure using AMPL.
- Flexible framework adaptable to various forest tree breeding scenarios.
- Separation of model and data for easy reuse and scalability.
- Case Study integration showing data generation and genetic evaluation in R.
- Compatible with multiple solvers (Gurobi is highly recommended).

---

## Repository Structure

### 1. AMPL Core Model (/AMPL)
Contains the fundamental algebraic modeling files.
- so.mod: Core model (parameters, variables, objective function, constraints).
- so.dat: Scalar parameters (e.g., number of candidates, status number).
- so.run: Execution script for loading the model, data, and running the optimization.

### 2. Case Study & R Scripts (/CaseStudy_R)
Demonstrates the practical application of the model using simulated data.
- 01_data_gen.R: Data generation using the MoBPS package and genetic parameter estimation using ASReml.
- 02_scenario_1.R: Implementation of optimization scenario 1.
- 03_scenario_2.R: Implementation of optimization scenario 2.

---

## Requirements
- AMPL: https://ampl.com
- Solver: A solver supporting MIQP (e.g., Gurobi).
- R Environment (Optional, for Case Study):
  - Required packages: MoBPS, ASReml-R

---

## How to Run

### Option A: Running the AMPL Model Directly
1. Install AMPL and a compatible solver.
2. Prepare your input data files (b.dat, l.dat, u.dat, c.dat) following the AMPL format.
3. Set scalar parameters in so.dat.
4. Execute the optimization in so.run

### Option B: Running the Case Study (R)
1. Run 01_data_gen.R to simulate the population and estimate genetic values (requires ASReml license).
2. Run 02_scenario_1.R or 03_scenario_2.R to execute the optimization across different scenarios. 

---

## Model Description
The optimization problem is defined as:
- Objective: Maximize genetic gain (weighted sum of breeding values and contributions).
- Constraints:
  - Total contribution equals 1.
  - Individual lower (l) and upper (u) bounds on contributions.
  - Constraint on Status Number (effective population size).
- Extensions: Supports pairwise coancestry, pollen contamination, and female/male contribution adjustments.

---

## Reproducibility
This repository provides all necessary components to reproduce the results described in the manuscript. Users can apply the model to their own datasets by modifying the input files in the AMPL/ folder or by adapting the R scripts provided in the CaseStudy_R/ directory.

---

### Contact
For any inquiries, please contact:

- Corresponding author: Prof. Milan Lstibůrek  
	E-mail: lstiburek@fld.czu.cz

- First author: Dr. Ye-Ji Kim  
	E-mail: kyeji1107@gmail.com
