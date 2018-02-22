################################################################################
# District.jl
################################################################################
# Implement PADP decomposition solver.
# --- Quantities decomposition scheme ---
################################################################################


export PADP


mutable struct PADP <: AbstractDecompositionSolver
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
    PADP(pb::Grid)

Build PADP solver.
"""
function PADP(pb::Grid; nsimu=100, nit=10, algo=SDDP(nit))
    nnodes = length(pb.nodes)
    ntime = ntimesteps(pb.nodes[1].time)

    λ = zeros(Float64, nnodes, ntime-1)
    μ = zeros(Float64, nnodes, ntime-1)
    scen = [genscen(d.model, nsimu) for d in pb.nodes]
    # initiate mod with empty dictionnary
    mod = Dict()

    PADP(ntime, Inf, λ, μ, algo, nothing, scen, nsimu, nit, mod)
end

################################################################################
# The three pilars of oracle
################################################################################

function solve!(pb::Grid, dadp::PADP)
    # solve production subproblems
    for d in pb.nodes
        dadp.models[d.name] = solve(d, dadp.algo)
    end
    # solve transport problem
    qsolve!(pb.net)
end

function simulate!(pb::Grid, dadp::PADP)
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

function ∇f(pb::Grid, dadp::PADP)
    dg = zeros(Float64, dadp.ntime-1, nnodes(pb))
    for t in 1:(dadp.ntime - 1)
        λ = dadp.λ[:, t]
        μ = dadp.μ[:, t]
        dg[t, :] = (λ - μ)[:]
    end
    return dg[:]
end
