################################################################################
# District.jl
################################################################################
# Implement ADMM decomposition solver.
# Design to solve graph problem, with interconnection balance.
#  min_{F, Q} J_P(F) + J_T(Q) ,
#      s.t. A Q + F = 0
################################################################################

export ADMM


struct ADMM <: AbstractSolver
    niter::Int
    approx::Symbol
    # stopping criterion
    rtol::Float64
    gtol::Float64
    nsimu::Int
    warmstart::Bool
    penalty::Float64
end

struct ADMMResults{T}
    minimizer::Array{T}
    minimum::T
    iterations::Int
    trace::Vector{T}
    f_calls::Int
    exectime::Vector{Float64}
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
        # reshape multiplier
        mul = reshape(λ, nnodes(pb), ntimes(pb)-1)

        swap!(pb, mul, F, Q)
        solve!(pb, solver)
        simulate!(pb, solver)

        # determine equilibrium along scenarios
        equilibrium = ∇g(pb)
        z = qflow(pb, Q)[:]

        # update multiplier
        λ += solver.penalty*equilibrium

        # stopping criterion
        ## primal residual
        ϵp = norm(equilibrium)
        ## dual residual
        ϵd = norm(solver.penalty*(z - zp) )

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

stopcrit(epsp, epsd, rtol, gtol)=(epsp < rtol) && (epsd < gtol)


# update Bellman values and estimate value with Monte-Carlo
function produpd!(pb::Grid, λ, Q, τ, nsimu)
    # TODO:
    nnodes, ntime, ninfo =  getsize(pb)
    f = zeros(Float64, nnodes, ntime-1, ninfo)
    for t in 1:ntime-1, ni in 1:ninfo
        q = Q[t, :, ni]
        f[:, t, ni] = pb.network.A*q
    end

    # TODO: update
    swap!(pb, λ, Q)
    # update Bellman functions
    for (i, node) in enumerate(pb.nodes)
        dpsolve!(node.problem, node.solver, node.V, atom,
                 node.solver.uspace, ADMMDecomposition())
    end

    # simulate
    F = zeros(Float64, nnodes, ntime-1, ninfo)
    for (i, node) in enumerate(pb.nodes)
        atom = ADMMAtom(λ[i, :, :], f[i, :, :,], τ)

        scenarios = StochDynamicProgramming.simulate_scenarios(node.problem.noises, nsimu)
        # TODO: dry simulation
        sim = simulate_uzawa(node.problem, node.solver, node.V,
                          node.x0, atom, node.solver.uspace, scenarios,
                          PriceDecomposition())
        F[i, :, :] = mean(sim[2], 2)
    end
    return F
end

function transupd!(pb::Grid, λ, F, τ)
    admmsolve!(pb.network, λ, F, τ)
    return pb.network.q
end
