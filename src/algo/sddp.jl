################################################################################
# District.jl
################################################################################
# Implement SDDP solver.
# - Currently, SDDP solver used StochDynamicProgramming.
################################################################################
# TODO: work in progress to interface efficiently StochDynamicProgramming

export SDDP

"""
    AbstractDPSolver

Abstract type to define Dynamic Programming solvers.
"""
abstract type AbstractDPSolver <: AbstractSolver end


# TODO: improve description of SDDP
immutable SDDP <: AbstractDPSolver
    nit::Int
end

# SOLVER
function solve(node::AbstractNode, sddp::SDDP)
    params = get_sddp_solver()
    return solve_SDDP(node.model, params, 0, 0,
                      stopcrit=IterLimit(sddp.nit))
end
