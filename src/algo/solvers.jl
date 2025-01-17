################################################################################
# District.jl
################################################################################
# Optimization solvers for Optimal Control.
# - Change here the LP solvers if you do not have Gurobi installed.
################################################################################
# TODO: clean SDDP solver

export solve

"""
    AbstractSolver

Abstract type to define optimization solver objects.
"""
abstract type AbstractSolver end


"""
    solve(model::AbstractNode, solver::AbstractSolver)

Solve problem specified in `model` with SP solver `solver`.
"""
function solve end


# UTILS
# Specify your solver:
using Gurobi
#= using Clp =#
# using CPLEX

#= get_solver() = Gurobi.GurobiSolver(OutputFlag=false, MIPGap=.01) =#
#get_solver() = CPLEX.CplexSolver(
#                                 CPX_PARAM_SIMDISPLAY=0,
#                                 CPX_PARAM_SCRIND=0,
#                                 CPX_PARAM_QPMETHOD=2
#                                )
"Get LP solver."
get_solver() = Gurobi.GurobiSolver(OutputFlag=false,
                                   Method=0,
                                   MIPGap=.01,
                                   Threads=1,
                                   TimeLimit=1
                                   )


"Get StochDynamicProgramming SDDP solver."
function get_sddp_solver()
    solver = get_solver()

    params = StochDynamicProgramming.SDDPparameters(solver,
                                                passnumber=1,
                                                montecarlo_in_iter=1000,
                                                montecarlo_final=-1,
                                                compute_ub=100,
                                                prune_cuts=-1,
                                                pruning_algo="exact+",
                                                gap=.001,
                                                max_iterations=300)
    params.monteCarloSize = 0
    return params
end
