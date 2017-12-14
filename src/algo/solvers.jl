
# Specify your solver:
using Gurobi
#= using Clp =#

"""Get solver."""
function get_solver()
    return Gurobi.GurobiSolver(OutputFlag=false, MIPGap=.01)
end


function get_sddp_solver(;mip=false)
    solver = get_solver()

    if mip
        solver2 =  Gurobi.GurobiSolver(OutputFlag=false, MIPGap=.01)
    else
        solver2 = nothing
    end

    params = StochDynamicProgramming.SDDPparameters(solver,
                                                passnumber=1,
                                                montecarlo_in_iter=1000,
                                                montecarlo_final=-1,
                                                compute_ub=100,
                                                prune_cuts=-1,
                                                mipsolver=solver2,
                                                pruning_algo="exact+",
                                                gap=.001,
                                                max_iterations=300)
    return params
end
