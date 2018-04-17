################################################################################
# District.jl
################################################################################
# Implement PADP decomposition solver.
# --- Quantities decomposition scheme ---
################################################################################


export QADP

# TODO: dry with PADP (almost similar)

mutable struct QADP <: AbstractDecompositionSolver
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
    QADP(pb::Grid)

Build QADP solver.
"""
function QADP(pb::Grid; nsimu=100, nit=10, algo=SDDP(nit))
    if ~checkconsistency(pb, FlowInterface)
        error("Wrong interfaces inside `pb.nodes`. Use `FlowInterface`
              for quantities decomposition")
    end
    nnodes = length(pb.nodes)
    na = size(pb.net.A, 2)
    ntime = ntimesteps(pb.nodes[1].time)

    λ = zeros(Float64, nnodes, ntime-1)
    μ = zeros(Float64, na, ntime-1)
    scen = [genscen(d.model, nsimu) for d in pb.nodes]
    # initiate mod with empty dictionnary
    mod = Dict()

    QADP(ntime, Inf, λ, μ, algo, scen, nsimu, nit, mod)
end


################################################################################
# The three pilars of oracle
################################################################################

function solve!(pb::Grid, dadp::QADP)
    # solve production subproblems
    for d in pb.nodes
        dadp.models[d.name] = solve(d, dadp.algo)
    end
    # solve transport problem
    flowsolve!(pb.net)
end

function simulate!(pb::Grid, dadp::QADP)
    # get sensitivity of production problem
    dadp.cost = dualsimulation!(pb, dadp)
    # add transportation cost
    dadp.cost += pb.net.cost
    # update Q flows inside DADP
    copy!(dadp.μ, pb.net.Q')
end

function ∇f(pb::Grid, dadp::QADP)
    dg = zeros(Float64, dadp.ntime-1, narcs(pb))
    for t in 1:(dadp.ntime - 1)
        λ = dadp.λ[:, t]
        μ = dadp.μ[:, t]
        # A' λ - μ
        dg[t, :] = (-pb.net.A' * λ + μ)[:]
    end
    return dg[:]
end
