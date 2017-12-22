
export SDDP

abstract type AbstractSolver end

# TODO: improve description of SDDP
immutable SDDP
    nit::Int
end

# SOLVER
function solve(node::AbstractNode, sddp::SDDP)
    params = District.get_sddp_solver()
    return solve_SDDP(node.model, params, 1,
                      stopcrit=District.IterLimit(sddp.nit))
end
