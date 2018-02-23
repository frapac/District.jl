################################################################################
# District.jl
################################################################################
# Implement policies to take online decisions in assessment.
# - MPCPolicy finds decisions with MPC algorithm.
# - DPPolicy finds decisions with a one-step lookahead problem.
################################################################################

export MPCPolicy

abstract type AbstractPolicy end

"""
    buildproblem!(p::AbstractPolicy, model::SPModel, t::Int)

Update JuMP.Model in `p.problem`, considering SPModel `model`
and current timestep `t`.
"""
function buildproblem! end

################################################################################
# MPC POLICY
################################################################################
abstract type AbstractMPCPolicy <: AbstractPolicy end

mutable struct MPCPolicy <: AbstractMPCPolicy
    # JuMP Model
    problem::JuMP.Model
    # Mathematical Programming solver
    solver
    # MPC forecast
    forecast
    # Horizon
    horizon
    # Final time
    final
end
MPCPolicy(forecast) = MPCPolicy(Model(), get_solver(), forecast, -1, size(forecast, 1))

# remaining time
ntimesteps(mpc::MPCPolicy, t) = mpc.final - t + 1

# TODO: clean
function build_oracle(forecast)
    function oracle(t)
        return collect(forecast[t, :])
    end
    return oracle
end


function buildproblem!(mpc::MPCPolicy, model, t::Int)
    oracle = build_oracle(mpc.forecast)
    ntime = ntimesteps(mpc, t)
    nx = model.dimStates
    nu = model.dimControls

    m = Model(solver=get_solver())

    @variable(m,  model.ulim[i][1] <= u[i=1:nu, j=1:ntime-1] <=  model.ulim[i][2])
    @variable(m,  model.xlim[i][1] <= x[i=1:nx, j=1:ntime] <= model.xlim[i][2])

    # TODO: clean final step
    @variable(m, zf >= 0)
    @constraint(m, zf >= 2. - x[2, ntime])

    @variable(m, w[1:model.dimNoises, 1:ntime])
    m.ext[:noise] = @constraint(m, w[:, 1] .== oracle(t))
    for j in 2:(ntime-1)
        @constraint(m, w[:, j] .== oracle(t+j-1))
    end

    cost0 = model.costFunctions(m, t, x[:, 1], u[:, 1], w[:, 1])
    costs = [model.costFunctions(m, t+j-1, x[:, j], u[:, j], w[:, j]) for j=2:ntime-1]

    # Set optimization objective:
    @objective(m, Min, cost0 + sum(costs[i] for i = 1:ntime-2) + PENAL_TANK*zf)

    @constraint(m, x[:, 2] .== model.dynamics(t, x[:, 1], u[:, 1], w[:, 1]))
    for j in 2:(ntime-1)
        # Dynamic constraint:
        @constraint(m, x[:, j+1] .== model.dynamics(j+t-1, x[:, j], u[:, j], w[:, j]))
    end

    # Add initial constraints:
    m.ext[:cons] = @constraint(m, init_cons, x[:, 1] .== 0)

    mpc.problem = m
end


# Take a decision with MPC
function (p::MPCPolicy)(x, ξ)
    m = p.problem
    u = m[:u]
    # update current state in MPC problem
    for i in 1:endof(x)
        JuMP.setRHS(m.ext[:cons][i], x[i])
    end
    # update current noise in MPC problem
    for i in 1:endof(ξ)
        JuMP.setRHS(m.ext[:noise][i], ξ[i])
    end
    # solve problem with LP solver
    st = JuMP.solve(m)
    # return first control
    return collect(getvalue(u)[:, 1])
end


"""Update forecast of MPC directly in JuMP Model."""
function updatemodel!(m, forecast, nscen)
    error("Deprecated")
    w = getvariable(m, :w)
    ntimes = size(w, 1)

    for k in 1:ntimes, m in 1:size(w, 2)
        JuMP.fix(w[m, k], forecast[m, nscen, k])
    end
end



################################################################################
# SDDP POLICY
################################################################################
abstract type AbstractDPPolicy <: AbstractPolicy end

ntime(p::AbstractDPPolicy) = length(p.V)

# Here and Now DP policy
mutable struct HereAndNowDP <: AbstractDPPolicy
    # JuMP Model
    problem::JuMP.Model
    # Value functions
    V
    # LP solver
    solver
end
HereAndNowDP(V) = HereAndNowDP(Model(), V, get_solver())

# Wait and See DP policy
# WaitAndSeeDP does not anticipate the realization of ξ_t between t and t+1,
# and thus take a decision wrt the law of ξ_t.
mutable struct WaitAndSeeDP <: AbstractDPPolicy
    # JuMP Model
    problem::JuMP.Model
    # Value functions
    V
    # LP solver
    solver
    # Laws of uncertainties {w_t}_{t=0..T}
    laws
end
WaitAndSeeDP(V, laws) = WaitAndSeeDP(Model(), V, get_solver(), laws)


function buildproblem!(policy::HereAndNowDP, model, t::Int)
    m = Model(solver=policy.solver)

    nx = model.dimStates
    nu = model.dimControls
    nw = model.dimNoises

    @variable(m,  model.xlim[i][1] <= x[i=1:nx] <= model.xlim[i][2])
    @variable(m,  model.xlim[i][1] <= xf[i=1:nx]<= model.xlim[i][2])
    @variable(m,  model.ulim[i][1] <= u[i=1:nu] <=  model.ulim[i][2],
              category=model.controlCat[i])
    @variable(m, alpha)

    @variable(m, w[1:nw] == 0)
    m.ext[:cons] = @constraint(m, state_constraint, x .== 0)

    @constraint(m, xf .== model.dynamics(t, x, u, w))

    @objective(m, Min, model.costFunctions(m, t, x, u, w) + alpha)

    for nc in 1:policy.V[t+1].numCuts
        lambda = vec(policy.V[t+1].lambdas[nc, :])
        @constraint(m, policy.V[t+1].betas[nc] + dot(lambda, xf) <= alpha)
    end

    if ~isnull(model.equalityConstraints)
        @constraint(m, get(model.equalityConstraints)(t, x, u, w) .== 0)
    end

    # Add final cost
    if t == ntime(policy) - 1
        model.finalCost(model, m)
    end

    policy.problem = m
end


function buildproblem!(policy::WaitAndSeeDP, model, t::Int)
    m = Model(solver=policy.solver)
    law = policy.laws

    nx = model.dimStates
    nu = model.dimControls
    nw = model.dimNoises

    ns = law[t].supportSize
    ξ = collect(law[t].support[:, :])
    πp = law[t].proba

    @variable(m, model.xlim[i][1] <= x[i=1:nx] <= model.xlim[i][2])
    @variable(m, model.ulim[i][1] <= u[i=1:nu] <=  model.ulim[i][2])
    @variable(m, model.xlim[i][1] <= xf[i=1:nx, j=1:ns]<= model.xlim[i][2])
    @variable(m, alpha[1:ns])

    @variable(m, w[1:nw] == 0)
    m.ext[:cons] = @constraint(m, state_constraint, x .== 0)

    for j=1:ns
        @constraint(m, xf[:, j] .== model.dynamics(t, x, u, ξ[:, j]))
    end

    # add objective as minimization of expectancy:
    @objective(m, Min,
               sum(πp[j]*(model.costFunctions(m, t, x, u, ξ[:, j]) +
                    alpha[j]) for j in 1:ns))

    for nc in 1:policy.V[t+1].numCuts
        lambda = vec(policy.V[t+1].lambdas[nc, :])
        for j=1:ns
            @constraint(m, policy.V[t+1].betas[nc] + dot(lambda, xf[:, j]) <= alpha[j])
        end
    end

    if ~isnull(model.equalityConstraints)
        @constraint(m, get(model.equalityConstraints)(t, x, u, w) .== 0)
    end

    policy.problem = m
end


function (p::HereAndNowDP)(x, ξ)
    m = p.problem
    u = m[:u]
    w = m[:w]

    # Update value of w:
    for ii in 1:endof(ξ)
        JuMP.fix(w[ii], ξ[ii])
    end
    # Update constraint x == xt
    for i in 1:endof(x)
        JuMP.setRHS(m.ext[:cons][i], x[i])
    end

    status = JuMP.solve(m, suppress_warnings=false)
    return getvalue(u)
end


# WaitAndSeeDP does not use the realization ξ to compute a decision
function (p::WaitAndSeeDP)(x, ξ)
    m = p.problem
    u = m[:u]

    # Update constraint x == xt
    for i in 1:endof(x)
        JuMP.setRHS(m.ext[:cons][i], x[i])
    end

    status = JuMP.solve(m, suppress_warnings=false)
    return getvalue(u)
end


################################################################################
# DADP POLICY
################################################################################
abstract type AbstractDADPPolicy <: AbstractPolicy end

ntime(p::AbstractDADPPolicy) = length(p.V[1])

# Here and Now DP policy
mutable struct DADPPolicy <: AbstractDADPPolicy
    # JuMP Model
    problem::JuMP.Model
    # Value functions
    V::Vector{Vector{PolyhedralFunction}}
    # LP solver
    solver
end
DADPPolicy(V) = DADPPolicy(Model(), V, get_solver())


function buildproblem!(policy::DADPPolicy, model, t::Int)
    m = Model(solver=policy.solver)

    # get number of value functions
    nvf = length(policy.V)
    nx = model.dimStates
    nu = model.dimControls
    nw = model.dimNoises

    @variable(m,  model.xlim[i][1] <= x[i=1:nx] <= model.xlim[i][2])
    @variable(m,  model.xlim[i][1] <= xf[i=1:nx]<= model.xlim[i][2])
    @variable(m,  model.ulim[i][1] <= u[i=1:nu] <=  model.ulim[i][2],
              category=model.controlCat[i])
    # TODO: get rid of alpha if not final Bellman operator
    @variable(m, alpha)

    @variable(m, w[1:nw] == 0)
    m.ext[:cons] = @constraint(m, state_constraint, x .== 0)

    @constraint(m, xf .== model.dynamics(t, x, u, w))

    # For all i ∈ [1, N], we set : α[i] = V_{t+1}^i (x_{t+1})
    @variable(m, α[1:nvf])
    # Here, we consider that global value function is
    # V^{tot}_t(x_t) = ∑_{i∈N} V_t^i(x_t^i)
    @objective(m, Min, model.costFunctions(m, t, x, u, w) + sum(α[i] for i in 1:nvf))

    xcount = 1
    for idim in 1:nvf
        Vi = policy.V[idim][t+1]
        nxi = size(Vi.lambdas, 2)

        for nc in 1:Vi.numCuts
            lambda = Vi.lambdas[nc, :]
            # we penalize only portion of xf corresponding to node `idim`
            @constraint(m, Vi.betas[nc] + dot(lambda, xf[xcount:xcount+nxi-1]) <= α[idim])
        end
        xcount += nxi
    end

    if ~isnull(model.equalityConstraints)
        @constraint(m, get(model.equalityConstraints)(t, x, u, w) .== 0)
    end

    # Add final cost
    if t == ntime(policy) - 1
        model.finalCost(model, m)
        # set nodal value functions to 0
        @constraint(m, α .== 0)
        @objective(m, Min, model.costFunctions(m, t, x, u, w) + alpha)
    end
    policy.problem = m
end


# DADPPolicy is basically a HereAndNow policy
function (p::DADPPolicy)(x, ξ)
    m = p.problem
    u = m[:u]
    w = m[:w]

    # Update value of w:
    for ii in 1:endof(ξ)
        JuMP.fix(w[ii], ξ[ii])
    end
    # Update constraint x == xt
    for i in 1:endof(x)
        JuMP.setRHS(m.ext[:cons][i], x[i])
    end

    status = JuMP.solve(m, suppress_warnings=false)
    return getvalue(u)
end
