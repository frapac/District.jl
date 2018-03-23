################################################################################
# District.jl
################################################################################
# Implement DADP decomposition solver.
# --- Price decomposition scheme ---
# Design to solve graph problem, with interconnection balance.
#  min_{F, Q} J_P(F) + J_T(Q) ,
#      s.t. A Q + F = 0
################################################################################

export DADP

import Base: Base.show

abstract type AbstractDecompositionSolver <: AbstractSolver end


# TODO: move params in dedicated structure
mutable struct DADP <: AbstractDecompositionSolver
    # number of timesteps
    ntime::Int
    # current dual cost
    cost::Float64
    # current flow in nodes
    F::Array{Float64, 2}
    # current flow in edges
    Q::Array{Float64}
    # solver to solve nodes subproblems
    algo::AbstractDPSolver
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
function DADP(pb::AbstractGrid; nsimu=100, nit=10, algo=SDDP(nit))
    if ~checkconsistency(pb, PriceInterface)
        error("Wrong interfaces inside `pb.nodes`. Use `PriceInterface`
              for price decomposition")
    end
    nnodes = nnodes(pb)
    ntime = ntimesteps(pb.nodes[1].time)

    F = zeros(Float64, nnodes, ntime-1)
    Q = zeros(Float64, ntime-1, pb.net.narcs)
    scen = [genscen(d.model, nsimu) for d in pb.nodes]
    # initiate mod with empty dictionnary
    mod = Dict()

    DADP(ntime, Inf, F, Q, algo, scen, nsimu, nit, mod)
end

# We need three ingredients to build oracle:
#     i)   solve! : once new multipliers are set, resolve nodes and edges
#                 problems and update value function.
#     ii)  simulate! : once value functions are obtained, simulate
#                     trajectories via Monte Carlo with relaxed problem.
#     iii) ∇g : eventually, compute subgradient along previous trajectories.

function solve!(pb::AbstractGrid, dadp::DADP)
    # solve production subproblems
    for d in pb.nodes
        dadp.models[d.name] = solve(d, dadp.algo)
    end
    # solve transport problem
    solve!(pb.net)
end

function simulate!(pb::Grid, dadp::DADP)
    dadp.cost = 0.
    for (id, d) in enumerate(pb.nodes)
        c, flow = mcsimulation(dadp.models[d.name], dadp.scen[id])
        # take average of importation flows for Node `d`
        dadp.F[id, :] = mean(flow, 2)
        # take average of costs
        dadp.cost -= mean(c)
    end

    # add transportation cost
    dadp.cost -= pb.net.cost
    # update Q flows inside DADP
    copy!(dadp.Q, pb.net.Q)
end

function simulate!(pb::ZonalGrid, dadp::DADP)
    dadp.cost = 0.
    nodeindex = 0.
    for (idzone, zone) in enumerate(pb.nodes)

        c, flow = mcsimulation(dadp.models[zone.name], dadp.scen[idzone], zone)
        for (id, d) in enumerate(zone.bordernodes)
            # take average of importation flows for Node `d`
            dadp.F[nodeindex + id, :] = mean(flow[id], 2)
        end

        # take average of costs
        dadp.cost -= mean(c)
        nodeindex += length(zone.bordernodes)
    end

    # add transportation cost
    dadp.cost -= pb.net.cost
    # update Q flows inside DADP
    copy!(dadp.Q, pb.net.Q)
end

function ∇f(pb::AbstractGrid, dadp::DADP)
    dg = zeros(Float64, dadp.ntime-1, nnodes(pb))
    for t in 1:(dadp.ntime - 1)
        f = dadp.F[:, t]
        q = dadp.Q[t, :]
        dg[t, :] = (pb.net.A*q + f)[:]
    end
    # minus sign because in DADP we consider -f instead of f  (max f = - min -f)
    return -dg[:]
end



################################################################################
# General oracle for decomposition
################################################################################
function _update!(pb::AbstractGrid, algo::AbstractDecompositionSolver, x::Vector{Float64})
    # update multiplier inside Grid `pb`
    swap!(pb, x)
    # resolve `pb` with these new multipliers
    solve!(pb, algo)
    # simulate trajectories with updated value functions
    simulate!(pb, algo)
end

# oracle return a cost function `f` and a gradient function `grad!`
# corresponding to the transporation problem.
function oracle(pb::AbstractGrid, algo::AbstractDecompositionSolver)
    xp = UInt64(0)
    function f(x)
        if hash(x) != xp
            _update!(pb, algo, x)
            xp = hash(x)
        end
        return algo.cost
    end

    function grad!(x, storage)
        if hash(x) != xp
            _update!(pb, algo, x)
            xp = hash(x)
        end
        # Then, compute subgradients!
        copy!(storage, ∇f(pb, algo))
    end

    return f, grad!
end
