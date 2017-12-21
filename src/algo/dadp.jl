
import Base: Base.show


    name::Symbol
    ntime::Int

    solver::SPSolver

    # multiplier
    λ::Array{Float64, 2}
    # value functions
    V::Array{Float64, 2}
    # optimal flow
    flow::Array{Float64, 2}

    # initial position
    x0::Vector{Float64}

    primalcost::Float64
    dualcost::Float64

    # size of information variable
    ninfo::Int


################################################################################
# Implementation of DADP algorithm
################################################################################


# TODO: add proper import
using Lbfgsb


type DADPAlgo
    N_ITER::Int
    algo::Optim.Optimizer
    approx::Symbol
    gtol::Float64
    nsimu::Int
    warmstart::Bool
end


"""DADP algorithm.

We use a solver to perform the gradient descent.
"""
function optimize!(pb::Grid, dadp::DADPAlgo; verbose=false)
    # get initial position
    x0 = init(pb.dec, pb, warmstart=dadp.warmstart)

    # build oracle
    f, grad! = oracle(pb, dadp.nsimu)

    # Configure Gradient Descent
    options = Optim.Options(g_tol=dadp.gtol, allow_f_increases=true, show_trace=verbose,
                        iterations=dadp.N_ITER, store_trace=true, extended_trace=true)

    # Launch Gradient Descent!
    try
        if isa(pb.dec, PriceDecomposition)
            #= return @time Optim.optimize(f, grad!, x0, dadp.algo, options) =#
            return @time lbfgsb(f, grad!, x0; iprint=1, pgtol=1e-8)
        end
    catch ex
        if isa(ex, InterruptException)
            warn("Gradient descent interrupted manually")
            return nothing
        else
            rethrow(ex)
        end
    end

end


# Define routine to use Gradient Descent.
# We use a closure to avoid unecessary operation
# while computing F after gradient G
function oracle(pb::Grid, nsimu::Int)
    function f(λ)
        return cost(pb)
    end

    function grad!(λ, storage)
        # update multiplier
        mul = reshape(λ, pb.nnodes, pb.ntime-1, pb.ninfo)

        isa(pb.dec, QuantDecomposition) && projection!(mul, pb)
        setexch!(pb, mul)

        # if decomposition is primal, use projected gradient
        isa(dec, QuantDecomposition) && proj!(pb, λ)

        # solve problem and update value functions
        solve!(pb)

        # Monte-Carlo simulation
        simulate!(pb, nsimu)

        # determine equilibrium along scenarios
        setsubgradient!(storage, pb)
    end

    return f, grad!
end


################################################################################
# SOLVER
################################################################################
function solve!(pb::Grid)
    # solve production problem
    psolve!(pb.dec, pb.nodes)
    # solve transport problem
    tsolve!(pb.dec, pb.network)
end


"""Solve transport problem."""
tsolve!(::PriceDecomposition, network::Network)=prsolve!(network)
tsolve!(::QuantDecomposition, network::Network)=qsolve!(network)

"""Solve production problem."""
function psolve!(dec::Decomposition, nodes::Vector{Node})
    #TODO: can be parallelized
    for node in nodes
        solve!(dec, node)
    end

end


################################################################################
# SIMULATION
################################################################################
function simulate!(pb::Grid, nsimu::Int)
    # TODO: parallelize simulation
    for node in pb.nodes
        seed = -1 #rand(1:1000)
        simnode!(pb.dec, node, nsimu, seed)
    end
end

function psimnode!(dec, node, nsimu)
    nt = Threads.nthreads()
    nsimu_th = div(nsimu, nt)

    cost = zeros(Float64, nt)
    F = zeros(Float64, size(node.flow)..., nt)
    Threads.@threads for i in 1:nt
        simd = simulate_uzawa(node, nsimu_th, PriceDecomposition(), seed)
        cost[i] = -mean(simd[1])
        F[:, :, i] = mean(simd[2], 2)
    end

    node.dualcost = mean(cost)
    # update flow inside node
    # TODO: mem alloc
    node.flow = mean(F, 3)[:, :, 1]
end

"""Simulate node in dual."""
function simnode!(::PriceDecomposition, node, nsimu, seed=-1)
    # perform simulation in dual
    simd = simulate_uzawa(node, nsimu, PriceDecomposition())
    # update dual cost
    node.dualcost = -mean(simd[1])
    # update flow inside node
    # TODO: mem alloc
    node.flow = mean(simd[2], 2)
end
"""Simulate node in primal."""
function simnode!(::QuantDecomposition, node, nsimu)
    # perform simulation in dual
    simd = simulate_uzawa(node, nsimu, QuantDecomposition())
    # update dual cost
    node.primalcost = mean(simd[1])
    # update flow inside node
    # TODO: mem alloc
    node.λ = mean(simd[2], 2)
end


"""Simulate global problem."""
function primalcost(pb::Grid, nsimu::Int=1000)
    return mean(simup(pb, nsimu)[1])
end
"""Monte-carlo simulation in primal."""
function simup(pb::Grid, nsimu::Int)
    # generate sample of scenarios
    al = StochDynamicProgramming.simulate_scenarios(pb.model.noises, nsimu)
    x0 = [n.x0[1] for n in pb.nodes]
    # simulate primal model
    return MPTS.simulate_primal(pb.model, pb.Δx, [n.V for n in pb.nodes], x0, al)
end


################################################################################
# SUBGRADIENT
################################################################################
# set multiple dispatch:
setsubgradient!(oracle, pb::Grid) = _setsubgradient!(oracle, pb, pb.dec)
_setsubgradient!(oracle, pb::Grid, ::PriceDecomposition) = equiflow!(oracle, pb)
_setsubgradient!(oracle, pb::Grid, ::QuantDecomposition) = equiprice!(oracle, pb)


################################################################################
# COST
################################################################################
cost(pb::Grid) = _cost(pb.dec, pb)
_cost(::PriceDecomposition, pb::Grid) = dualcost(pb)
_cost(::QuantDecomposition, pb::Grid) = primcost(pb)


"""Get evolution of primal cost through iterations."""
function getprimal(pb, trace, nsimu)
    primal = []
    dec = MPTS.PriceDecomposition()

    for t in 1:endof(trace)
        _mul = trace[t].metadata["x"]
        mul = reshape(_mul, pb.nnodes, pb.ntime-1, pb.ninfo)
        MPTS.setexch!(pb, mul)

        MPTS.psolve!(dec, pb.nodes)
        prim = MPTS.primalcost(pb, nsimu)
        push!(primal, prim)
    end
    primal
end

