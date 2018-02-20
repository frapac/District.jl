################################################################################
# District.jl
################################################################################
# Solve transport problems.
################################################################################
# TODO: clean different solve! function with type dispatch

"""Solve transport problem in dual."""
function solve!(net::Network)
    tcost = 0.
    k1 = net.k1
    k2 = net.k2
    λ = net.λ
    flows = zeros(Float64, net.ntime-1, size(net.A, 2))

    for t in 1:net.ntime-1

        # TODO: dry definition of problem to use hotstart
        m = Model(solver=get_solver())
        mul = @view λ[t, :]
        @variable(m, -net.maxflow[i] <= q[i=1:net.narcs] <= net.maxflow[i])

        # take qp = |q|
        @variable(m, qp[1:net.narcs])
        @constraint(m, qp .>=  q)
        @constraint(m, qp .>= -q)

        @objective(m, :Min, sum(k1*qp)  + k2*dot(q, q) +  dot(mul, net.A*q))

        # solve problem
        JuMP.solve(m)
        qa = getvalue(q)
        tcost += getobjectivevalue(m)
        flows[t, :] = qa
    end

    net.cost = tcost
    net.Q = flows
end

"""Solve transport problem in primal."""
function qsolve!(net::Network, flow)
    k1, k2 = net.k1, net.k2

    narcs = net.narcs
    qa = zeros(Float64, net.ntime-1, size(net.A, 2))
    pcost = 0.

    for t in 1:net.ntime-1
        f = @view flow[:, t]
        # TODO: dry definition of problem to use hotstart
        m = Model(solver=get_solver())
        @variable(m, -net.maxflow[i] <= q[i=1:narcs] <= net.maxflow[i])

        # equilibrium
        eq = @constraint(m, net.A*q + f .== 0)

        # take qp = |q|
        @variable(m, qp)
        @constraint(m, qp .>=  q)
        @constraint(m, qp .>= -q)

        @objective(m, :Min, k1*sum(qp) + k2*dot(q, q))

        # solve problem
        status = JuMP.solve(m)
        if status != :Optimal
            println(sum(f))
            println(m)
        end
        @assert status == :Optimal

        qa[t, :] = getvalue(q)

        # update multiplier
        net.λ[t, :] = getdual(eq)
        pcost += getobjectivevalue(m)
    end

    # update q
    net.cost = pcost
    net.Q = qa
end


function admmsolve!(net::Network, λ, F, τ)
    k1, k2 = net.k1, net.k2

    narcs = net.narcs
    qa = zeros(Float64, net.ntime-1, size(net.A, 2))
    pcost = 0.

    for t in 1:net.ntime-1
        mul = @view λ[:, t]
        f =   @view F[:, t]

        m = Model(solver=get_solver())
        @variable(m, -net.maxflow[i] <= q[i=1:narcs] <= net.maxflow[i])

        # take qp = |q|
        @variable(m, qp)
        @constraint(m, qp .>=  q)
        @constraint(m, qp .>= -q)

        @objective(m, :Min, sum(k1*qp) + k2*dot(q, q) +
                            dot(mul, net.A*q + f) + τ*dot(net.A*q+f, net.A*q+f))

        # solve problem
        status = JuMP.solve(m)

        qa[t, :] = getvalue(q)

        # update multiplier
        pcost += getobjectivevalue(m)
    end

    net.cost = pcost
    # update q
    net.Q = qa
end
