# District.jl

A stochastic optimization library for district microgrids.


District.jl is a modeling and optimization library, aiming at
applying stochastic optimization algorithms to the control
of local microgrids. The three components of District are:

- A modeling tool to build easily optimization model;
- A set of optimization algorithms to solve these optimization models;
- A simulator to fairly assess the different strategies.

More information in documentation.


## Nota-Bene

Current version works with:

- Julia 0.6
- JuMP 0.18.2 (f0bf187)
- Gurobi 0.4.1 (f54387f)
- StochDynamicProgramming 0.5.0 (d9b339d)

Integration with Julia 1.0 remains to be done.
