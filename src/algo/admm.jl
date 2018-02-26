################################################################################
# District.jl
################################################################################
# Implement ADMM decomposition solver.
# Design to solve graph problem, with interconnection balance.
#  min_{F, Q} J_P(F) + J_T(Q) ,
#      s.t. A Q + F = 0
################################################################################

export ADMM


struct ADMM <: AbstractDecompositionSolver
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
    # quad penalty term
    τ::Float64
    # stopping criterion
    rtol::Float64
    gtol::Float64
    nsimu::Int
    nit::Int
    models::Dict()
end
"""
    ADMM(pb::Grid)

Build ADMM solver.
"""
function ADMM(pb::Grid; nsimu=100, nit=10, algo=SDDP(nit), rtol=1e-2, gtol=1e-2, tau=1.)
    if ~checkconsistency(pb, QuadInterface)
        error("Wrong interfaces inside `pb.nodes`. Use `PriceInterface`
              for price decomposition")
    end
    nnodes = length(pb.nodes)
    ntime = ntimesteps(pb.nodes[1].time)

    F = zeros(Float64, nnodes, ntime-1)
    Q = zeros(Float64, ntime-1, pb.net.narcs)
    scen = [genscen(d.model, nsimu) for d in pb.nodes]
    # initiate mod with empty dictionnary
    mod = Dict()

    ADMM(ntime, Inf, F, Q, algo, scen, tau, rtol, gtol, nsimu, nit, mod)
end

struct ADMMResults{T}
    minimizer::Array{T}
    minimum::T
    iterations::Int
    trace::Vector{T}
    f_calls::Int
    exectime::Vector{Float64}
end

function solve!(pb::Grid, dadp::ADMM)
    # solve production subproblems
    for d in pb.nodes
        dadp.models[d.name] = solve(d, dadp.algo)
    end
    # solve transport problem with augmented Lagrangian
    admmsolve!(pb.net, dadp.F, dadp.τ)
end

function simulate!(pb::Grid, dadp::ADMM)
    dadp.cost = 0.
    for (id, d) in enumerate(pb.nodes)
        c, _, u = StochDynamicProgramming.simulate(dadp.models[d.name], dadp.scen[id])
        # take average of importation flows for Node `d`
        dadp.F[id, :] = mean(u[:, :, end], 2)
        # take average of costs
        dadp.cost -= mean(c)
    end

    # add transportation cost
    dadp.cost -= pb.net.cost
    # update Q flows inside DADP
    copy!(dadp.Q, pb.net.Q)
end

function ∇f(pb::Grid, dadp::ADMM)
    dg = zeros(Float64, dadp.ntime-1, nnodes(pb))
    for t in 1:(dadp.ntime - 1)
        f = dadp.F[:, t]
        q = dadp.Q[t, :]
        dg[t, :] = (pb.net.A*q + f)[:]
    end
    # minus sign because in DADP we consider -f instead of f  (max f = - min -f)
    return -dg[:]
end



"""ADMM algorithm.

We use a solver to perform the gradient descent.
"""
function solve!(pb::Grid, xini::Vector{Float64}, solver::ADMM; verbose=false)
    # get initial position
    λ = copy(xini)
    equilibrium = zeros(Float64, length(λ))

    Q = copy(pb.net.Q)
    λp = Inf
    Qp = Inf
    Fp = Inf
    zp = Inf

    # store evolution of gradient
    it = -1
    ngrads = Float64[]
    exec = Float64[]

    for it in 1:solver.niter
        tic()
        # TODO: add multiswap for grids
        swap!(pb, λ, F)
        solve!(pb, solver)
        simulate!(pb, solver)

        # determine equilibrium along scenarios
        equilibrium = ∇g(pb)
        z = qflow(pb, Q)[:]

        # update multiplier
        λ -= solver.τ*equilibrium

        # stopping criterion
        ## primal residual
        ϵp = norm(equilibrium)
        ## dual residual
        ϵd = norm(solver.τ *(z - zp))

        stopcrit(ϵp, ϵd, solver.rtol, solver.gtol) && break

        tf = toq()
        verbose && @printf("\t %s \t  %.4e \t %.4e \t %.4e \n", it, ϵp, ϵd, tf)

        push!(ngrads, norm(equilibrium))
        push!(exec, tf)

        # save current values
        λp = copy(λ)
        Fp = copy(F[:])
        Qp = copy(Q[:])
        zp = copy(z)
    end

    ADMMResults(λ, -Inf, it, ngrads, it, exec)
end

stopcrit(epsp, epsd, rtol, gtol) = (epsp < rtol) && (epsd < gtol)
