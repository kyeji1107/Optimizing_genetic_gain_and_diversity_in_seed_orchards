### Optimization Framework using AMPL
This repository provides the AMPL implementation of the optimization framework presented in the manuscript:
“Optimizing genetic gain and diversity in seed orchards using an algebraic modeling language for mathematical programming.”
The framework formulates seed orchard deployment as a mixed-integer quadratic programming (MIQP) problem to maximize genetic gain while maintaining constraints on genetic diversity.

### Overview
- Algebraic and transparent model structure using AMPL
- Flexible framework adaptable to different breeding scenarios
- Separation of model and data for easy reuse
- Compatible with multiple solvers (Gurobi recommended)

### Repository Structure (need to be revised)
- so.mod  
	Core AMPL model (parameters, variables, objective function, constraints)
- so.dat  
	Scalar parameters (e.g., number of candidates, status number)
- b.dat  
	Breeding values (selection criterion)
- l.dat, u.dat  
	Lower and upper bounds for contributions
- c.dat  
	Coancestry or genomic relationship matrix
- mobps.r
  	Execution script #1 for the case studies - generating simulated populations 
- so.r  
	Execution script #2 for the case studies - building model, and running optimization

### Model Description
The optimization problem is defined as:
- Objective:
	- Maximize genetic gain (weighted sum of breeding values and contributions)
- Constraints:
	- Total contribution equals 1
	- Bounds on individual contributions
	- Constraint on Status Number
- Extensions (optional):
	- Pairwise coancestry constraints
	- Pollen contamination
	- Female and male contribution

### Requirements
- AMPL (https://ampl.com)
- Solver supporting MIQP:
	- Gurobi (recommended)
	- Other compatible solvers may also be used

### Data Preparation
- Input data must follow the format expected by AMPL
- Data can be generated using R, Python, MATLAB, or other tools
- The framework is independent of the data generation process

### Reproducibility
This repository provides all necessary components to reproduce the optimization framework described in the manuscript. Users can apply the model to their own datasets by modifying input files without changing the core model.

### Contact
For any inquiries, please contact:

- Corresponding author: Prof. Milan Lstibůrek  
	Email: lstiburek@fld.czu.cz

- First author: Dr. Ye-Ji Kim  
	Email: kyeji1107@gmail.com
