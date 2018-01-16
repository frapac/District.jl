################################################################################
# District.jl
################################################################################
# Implement DADP decomposition solver.
# Design to solve graph problem, with interconnection balance.
#  min_{F, Q} J_P(F) + J_T(Q) ,
#      s.t. A Q + F = 0
################################################################################

export DADP

import Base: Base.show



# TODO: move params in dedicated structure
mutable struct DADP <: AbstractSolver
    # number of timesteps
    ntime::Int
    # current dual cost
    cost::Float64
    # current flow in nodes
    F::Array{Float64, 2}
    # current flow in edges
    Q::Array{Float64}
    # solver to solve nodes subproblems
    algo::AbstractSolver
    # ???
    pb
    # Scenario
    scen::Array
    # number of Monte Carlo simulations to estimate gradient
    nsimu::Int
    # maximum number of iterations
    nit::Int
    # models stores SDDPInterface in dedicated Dictionnary
    models::Dict
end


"""
    DADP(pb::Grid)

Build DADP solver.
"""
function DADP(pb::Grid; nsimu=100, nit=10, algo=SDDP(nit))
    nnodes = length(pb.nodes)
    ntime = ntimesteps(pb.nodes[1].time)

    F = zeros(Float64, nnodes, ntime-1)
    Q = zeros(Float64, ntime-1, pb.net.narcs)
    scen = [genscen(d.model, nsimu) for d in pb.nodes]
    # initiate mod with empty dictionnary
    mod = Dict()

    DADP(ntime, Inf, F, Q, algo, nothing, scen, nsimu, nit, mod)
end

# We need three ingredients to build oracle:
#     i) solve! : once new multipliers are set, resolve nodes and edges
#                 problems and update value function.
#     ii) simulate! : once value functions are obtained, simulate
#                     trajectories via Monte Carlo with relaxed problem.
#     iii) ∇g : eventually, compute subgradient along previous trajectories.

function solve!(pb::Grid, dadp::DADP)
    # solve production subproblems
    for d in pb.nodes
        dadp.models[d.name] = solve(d, dadp.algo)
    end
    # solve transport problem
    solve!(pb.net)
end

function simulate!(pb::Grid, dadp::DADP)
    dadp.cost = 0
    for (id, d) in enumerate(pb.nodes)
        c, _, u = StochDynamicProgramming.simulate(dadp.models[d.name], dadp.scen[id])
        # take average of importation flows for Node `d`
        dadp.F[id, :] = mean(u[:, :, end], 2)
        # take average of costs
        dadp.cost += mean(c)
    end

    # add transportation cost
    dadp.cost += pb.net.cost
    # update Q flows inside DADP
    copy!(dadp.Q, pb.net.Q)
end

function ∇g(pb::Grid, dadp::DADP)
    dg = zeros(Float64, dadp.ntime-1, nnodes(pb))
    for t in 1:(dadp.ntime - 1)
        f = dadp.F[:, t]
        q = dadp.Q[t, :]
        dg[t, :] = (pb.net.A*q + f)[:]
    end
    return dg[:]
end

# oracle return a cost function `f` and a gradient function `grad!`
# corresponding to the transporation problem.
function oracle(pb::Grid, dadp::DADP)
    # take care: we aim at find a maximum (min f = - max -f )
    f(λ) = - dadp.cost

    function grad!(λ, storage)
        # update multiplier inside Grid `pb`
        swap!(pb, λ)
        # resolve `pb` with these new multipliers
        solve!(pb, dadp)
        # simulate trajectories with updated value functions
        simulate!(pb, dadp)

        # Then, compute subgradients!
        copy!(storage, -∇g(pb, dadp))
    end

    return f, grad!
end
