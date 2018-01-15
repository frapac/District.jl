################################################################################
# District.jl
################################################################################
# Implement SDDP solver.
# - Currently, SDDP solver used StochDynamicProgramming.
################################################################################
# TODO: work in progress to interface efficiently StochDynamicProgramming

export SDDP


# TODO: improve description of SDDP
immutable SDDP <: AbstractSolver
    nit::Int
end

# SOLVER
function solve(node::AbstractNode, sddp::SDDP)
    params = District.get_sddp_solver()
    return solve_SDDP(node.model, params, 1,
                      stopcrit=District.IterLimit(sddp.nit))
end
