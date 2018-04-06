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
    V::Array{Float64}
    # solver to solve nodes subproblems
    algo::AbstractDPSolver
    # Scenario
    scen::Array
    # number of Monte Carlo simulations to estimate gradient
    nsimu::Int
    # maximum number of iterations
    maxit::Int
    # models stores SDDPInterface in dedicated Dictionnary
    models::Dict
end


"""
    MADP(pb::Grid)

Build MADP solver.
"""
function MADP(pb::Grid; nsimu=100, nit=10, algo=SDDP(nit), maxit=50)
    if ~checkconsistency(pb, FlowInterface)
        error("Wrong interfaces inside `pb.nodes`. Use `FlowInterface`
              for prediction decomposition")
    end
    nnodes = length(pb.nodes)
    ntime = ntimesteps(pb.nodes[1].time)

    λ = zeros(Float64, nnodes, ntime-1)
    μ = zeros(Float64, nnodes, ntime-1)
    scen = [genscen(d.model, nsimu) for d in pb.nodes]
    # initiate mod with empty dictionnary
    mod = Dict()

    MADP(ntime, Inf, λ, μ, algo, scen, nsimu, maxit, mod)
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

function simulate!(pb::Grid, dadp::MADP)
    # production cost
    dadp.cost = dualsimulation!(pb, dadp)
    # add transportation cost
    dadp.cost += pb.net.cost
end

flowallocation(pb::Grid) = flowallocation(pb.net)

# Solve problem with a fixed point algorithm.
function solve!(pb::Grid, algo::MADP, V0::Vector{Float64}, μ0::Vector{Float64})
    # alloc
    V = copy(V0)
    # price
    μ = copy(μ0)

    for nit = 1:algo.maxit
        prodswap!(pb, V)
        transswap!(pb, μ)
        # resolve
        solve!(pb, algo)
        # simulate trajectories with updated value functions
        simulate!(pb, algo)

        # update allocation and transport price
        V[:] = -flowallocation(pb)
        println(mean(pb.nodes[1].conn.values))
        μ[:] = (algo.λ')[:]

        @printf("\t %i \t %.6e\n", nit, algo.cost)
    end
    return V, μ
end
