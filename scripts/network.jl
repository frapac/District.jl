
mutable struct Network
    ntime::Int
    cost::Float64
    Q::Array{Float64, 2}
    # incidence matrix
    A::Array{Float64, 2}
    narcs::Int
    maxflow::Vector{Float64}

    # multiplier
    λ::Array{Float64, 2}
end
function Network(ts, A)
    ntime = District.ntimesteps(ts)
    nnodes, narcs = size(A)
    λ = zeros(Float64, ntime-1, nnodes)
    Q = zeros(Float64, ntime-1, narcs)
    # TODO: dry
    maxflow = 6 * ones(Float64, narcs)

    return Network(ntime, 1e3, Q, A, narcs, maxflow, λ)
end

swap!(net::Network, mul) = net.λ[:] = mul


"""Solve transport problem in dual."""
function solve!(net::Network)
    tcost = 0.
    k1 = 0.01
    #= λ = reshape(net.λ, net.ntime-1, size(net.A, 1), net.ninfo) =#
    λ = net.λ
    flows = zeros(Float64, net.ntime-1, size(net.A, 2))

    for t in 1:net.ntime-1

        # TODO: dry definition of problem to use hotstart
        m = Model(solver=GurobiSolver(OutputFlag=0))
        mul = @view λ[t, :]
        @variable(m, -net.maxflow[i] <= q[i=1:net.narcs] <= net.maxflow[i])

        # take qp = |q|
        @variable(m, qp)
        @constraint(m, qp .>=  q)
        @constraint(m, qp .>= -q)

        @objective(m, :Min, sum(k1*qp)  + dot(q, q) +  dot(mul, net.A*q))

        # solve problem
        JuMP.solve(m)
        qa = getvalue(q)
        tcost += getobjectivevalue(m)
        flows[t, :] = qa
    end

    net.cost = tcost
    net.Q = flows
end
