# Documentation of District.jl

District.jl is a package implementing various algorithms for
the controls of microgrid.


Please find attached a first version of the documentation
of District.jl package.


District.jl structure is divided in three parts:
- First, `src/model/` implements the modeling of the different devices we consider in local microgrid, and metaprogramming tools to build optimization model on the fly.
- Second, `src/algo` implements different algorithms to solve the optimization problems.
- Eventually, `src/simulation` implements simulation tools to fairly assess the different algorithms.

```@meta
CurrentModule = District
```
```@contents
Pages = ["model.md", "solver.md", "simulation.md"]
Depth = 3
```
