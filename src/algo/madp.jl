################################################################################
# District.jl
################################################################################
# Implement MADP decomposition solver.
# --- Decomposition by prediction scheme ---
################################################################################


export MADP


mutable struct MADP <: AbstractDecompositionSolver
    # number of timesteps
    ntime::Int
    # current dual cost
    cost::Float64
    # current price / nodes
    λ::Array{Float64, 2}
    # current price / edges
    μ::Array{Float64}
    # solver to solve nodes subproblems
    algo::AbstractDPSolver
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
    MADP(pb::Grid)

Build MADP solver.
"""
function MADP(pb::Grid; nsimu=100, nit=10, algo=SDDP(nit))
    nnodes = length(pb.nodes)
    ntime = ntimesteps(pb.nodes[1].time)

    λ = zeros(Float64, nnodes, ntime-1)
    μ = zeros(Float64, nnodes, ntime-1)
    scen = [genscen(d.model, nsimu) for d in pb.nodes]
    # initiate mod with empty dictionnary
    mod = Dict()

    MADP(ntime, Inf, λ, μ, algo, nothing, scen, nsimu, nit, mod)
end

################################################################################
# The three pilars of oracle
################################################################################

function solve!(pb::Grid, dadp::MADP)
    # solve production subproblems
    for d in pb.nodes
        dadp.models[d.name] = solve(d, dadp.algo)
    end
    # solve transport problem in dual
    solve!(pb.net)
end

# TODO: dry
function simulate!(pb::Grid, dadp::MADP)
    dadp.cost = 0.
    for (id, d) in enumerate(pb.nodes)
        c, λ = qsensitivity(dadp.models[d.name], dadp.scen[id])
        # take average of sensitivity wrt flows for Node `d`
        avgλ = mean(λ[:, :, end], 2)
        # sometimes, inner solver is unable to return the dual value
        # corresponding to the coupling constraint, and return NaN.
        # When this happened, we replace the NaN values by their
        # previous estimation, stored in dadp.λ
        avgλ[isnan.(avgλ)] = dadp.λ[id, vec(isnan.(avgλ))]
        # store new sensitivity
        dadp.λ[id, :] = avgλ
        # take average of costs
        dadp.cost += mean(c)
    end

    # add transportation cost
    dadp.cost += pb.net.cost
    # update Q flows inside DADP
    copy!(dadp.μ, pb.net.F')
end


# Solve problem with a fixed point algorithm.
function solve!(pb::Grid, algo:MADP, x0::Vector{Float64})

    # alloc
    V
    # price

    for ???
        swap!(pb, x)
        # resolve
        solve!(pb, algo)
        # simulate trajectories with updated value functions
        simulate!(pb, algo)

        V = A*pb.net.Q

    end

    return res
end

