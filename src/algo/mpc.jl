
"""Solve linear problem at each iteration of MPC.

# Arguments
* `mpc::MPC`
* `t::Int`
    Current time
* `oracle`
    Forecast to consider
* `x0`
    State of the system at time `t`
* `Tf::Int`
    Final time
"""
function mpcproblem(model, mpc, t, oracle, Tf)

    if mpc.SHRINK
        ntime = Tf - t + 1
    else
        ntime = mpc.horizon
    end
    nx = model.dimStates
    nu = model.dimControls

    m = Model(solver=get_solver())

    # take into account binary constraints in MPC:
    @variable(m,  model.ulim[i][1] <= u[i=1:nu, j=1:ntime-1] <=  model.ulim[i][2])
    # Define constraints other variables:
    #TODO: fix final time
    @variable(m,  model.xlim[i][1] <= x[i=1:nx, j=1:ntime] <= model.xlim[i][2])

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

    return m
end




function mpccontrol(model, m, x0)
    u = m[:u]
    for i in 1:model.dimStates
        JuMP.setRHS(m.ext[:cons][i], x0[i])
    end
    st = solve(m)
    if st == :Optimal
        return collect(getvalue(u)[:, 1])
    else
        return zeros(model.dimControls)
    end
end


function mpccontrol(model, m, x0, w0)
    for i in 1:model.dimNoises
        JuMP.setRHS(m.ext[:noise][i], w0[i])
    end
    try
        return mpccontrol(model, m, x0)
    catch
        println( w0)
    end
end

"""Update forecast of MPC directly in JuMP Model."""
function updatemodel!(m, forecast, nscen)
    w = getvariable(m, :w)
    ntimes = size(w, 1)

    for k in 1:ntimes, m in 1:size(w, 2)
        JuMP.fix(w[m, k], forecast[m, nscen, k])
    end
end
