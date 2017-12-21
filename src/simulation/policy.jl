################################################################################
# MPC POLICY
################################################################################

abstract type AbstractMPCPolicy end

struct MPCPolicy <: AbstractMPCPolicy
    problem
    solver
    forecast
end

# TODO: clean
function build_oracle(forecast)
    function oracle(t)
        return collect(forecast[t, :])
    end
    return oracle
end


"""Solve linear problem at each iteration of MPC.  """
function buildproblem!(mpc::MPCPolicy, model, t::Int)
    oracle = build_oracle(mpc.forecast)
    ntime = ntimesteps(mpc)
    nx = model.dimStates
    nu = model.dimControls

    m = Model(solver=get_solver())

    # take into account binary constraints in MPC:
    @variable(m,  model.ulim[i][1] <= u[i=1:nu, j=1:ntime-1] <=  model.ulim[i][2])
    # Define constraints other variables:
    @variable(m,  model.xlim[i][1] <= x[i=1:nx, j=1:ntime] <= model.xlim[i][2])

    # TODO: clean final step
    @variable(m, zf >= 0)
    @constraint(m, zf >= 6 - x[2, ntime])

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


function (p::MPCPolicy)(x, ξ)
    m = p.problem
    u = m[:u]
    for i in 1:model.dimStates
        JuMP.setRHS(m.ext[:cons][i], x[i])
    end
    for i in 1:model.dimNoises
        JuMP.setRHS(m.ext[:noise][i], ξ[i])
    end
    st = solve(m)
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
abstract type AbstractDPPolicy end

struct HereAndNowDP <: AbstractDPPolicy
    problem::JuMP.Model
    V
    solver
end
struct WaitAndSeeDP <: AbstractDPPolicy
    problem::JuMP.Model
    V
    solver
    laws
end


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

    if isa(model.costFunctions, Function)
    try
        @objective(m, Min, model.costFunctions(t, x, u, w) + alpha)
    catch
        @objective(m, Min, model.costFunctions(m, t, x, u, w) + alpha)
    end

    for nc in 1:policy.V[t+1].numCuts
        lambda = vec(V[t+1].lambdas[nc, :])
        @constraint(m, V[t+1].betas[nc] + dot(lambda, xf) <= alpha)
    end
    # TODO: fix final costs in definition of costs
    policy.problem = m
end


function buildproblem!(policy::WaitAndSeeDP, model, t)
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

    for nc in 1:V[t+1].numCuts
        lambda = vec(V[t+1].lambdas[nc, :])
        for j=1:ns
            @constraint(m, V[t+1].betas[nc] + dot(lambda, xf[:, j]) <= alpha[j])
        end
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

    status = solve(m, suppress_warnings=false)
    return getvalue(u)
end


function (p::WaitAndSeeDP)(x, ξ)
    m = p.problem
    u = m[:u]

    # Update constraint x == xt
    for i in 1:endof(x)
        JuMP.setRHS(m.ext[:cons][i], x[i])
    end

    status = solve(m, suppress_warnings=false)
    return getvalue(u)
end
