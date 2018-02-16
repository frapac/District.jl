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
