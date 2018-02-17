################################################################################
# District.jl
################################################################################
# Solve transport problems.
################################################################################

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
function qsolve!(net::Network)
    k1, k2 = net.k1, net.k2

    narcs = net.narcs
    qa = zeros(Float64, net.ntime-1, size(net.A, 2))
    pcost = 0.

    for t in 1:net.ntime-1, ni in 1:net.ninfo
        f = @view net.flow[:, t, ni]
        # TODO: dry definition of problem to use hotstart
        m = Model(solver=getsolver())
        @variable(m, -net.maxflow[i] <= q[i=1:narcs] <= net.maxflow[i])

        # equilibrium
        eq = @constraint(m, net.A*q + f .== 0)

        # take qp = |q|
        @variable(m, qp)
        @constraint(m, qp .>=  q)
        @constraint(m, qp .>= -q)

        @objective(m, :Min, sum(k1*qp + k2*dot(q, q)))

        # solve problem
        status = JuMP.solve(m)
        if status != :Optimal
            println(sum(f))
            println(m)
        end
        @assert status == :Optimal

        qa[t, :, ni] = getvalue(q)

        # update multiplier
        net.λ[:, t, ni] = getdual(eq)
        pcost += getobjectivevalue(m)
    end

    # update q
    net.cost = pcost
    net.Q = qa
end
