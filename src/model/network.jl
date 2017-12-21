

type Network
    ntime::Int
    # incidence matrix
    A::Array{Float64, 2}
    narcs::Int
    maxflow::Vector{Float64}

    # dual cost
    dualcost::Float64
    # primal cost
    primalcost::Float64

    # multiplier
    λ::Array{Float64, 3}

    # flow through arcs
    # size = (ntime, narcs, ny)
    q::Array{Float64, 3}
    flow::Array{Float64, 3}


    ninfo::Int
end
function Network(ntime, names::Vector{Symbol}, ny::Int)
    pos = getpos.(names)

    nstages = (ntime>0) ? ntime : Configuration.NTIME
    κ = Configuration.NTIME / nstages

    # we restrict connexion to selected position
    connex = build_graph()[pos, pos]

    # get number of nodes
    nnodes = size(connex, 1)
    @assert nnodes == length(names)
    # get number of arcs
    narcs = floor(Int, sum(connex .> 0.)/2)

    A, mf = buildincidence(connex)

    # For the sake of homogeneity, we write multipliers as vector
    λ = zeros(Float64, nnodes,(ntime-1),ny)
    # and flows as tensor
    q = zeros(Float64, ntime-1, narcs, ny)
    flow = zeros(Float64, ntime-1, nnodes, ny)
    Network(ntime, A, narcs, κ*mf, Inf, Inf, λ, q, flow, ny)
end


"""Solve network problem and return admissible flow `Q`."""
function solve!(net::Network)

    # we compute analytically the optimal flow
    k1, k2 = MPTS.Configuration.κ1, MPTS.Configuration.κ2
    k2 *= k1

    tcost = 0.
    #= λ = reshape(net.λ, net.ntime-1, size(net.A, 1), net.ninfo) =#
    λ = net.λ
    q = zeros(Float64, net.ntime-1, size(net.A, 2), net.ninfo)

    for t in 1:net.ntime-1, ni in 1:net.ninfo
        mul = @view λ[:, t, ni]
        qa = 1/(2*k2)*(-net.A'*mul - k1)
        tcost += sum(k1*qa + k2*qa.^2 )+ mul'*net.A*qa
        q[t, :, ni] = qa
    end

    #TODO: mem alloc
    net.q = q
    net.dualcost = -tcost[1]
end


function admmsolve!(net::Network, λ, F, τ)
    k1, k2 = MPTS.Configuration.κ1, MPTS.Configuration.κ2
    k2 *= k1

    narcs = net.narcs
    qa = zeros(Float64, net.ntime-1, size(net.A, 2), net.ninfo)
    pcost = 0.::Float64

    for t in 1:net.ntime-1, ni in 1:net.ninfo
        mul = @view λ[:, t, ni]
        f =   @view F[:, t, ni]
        m = Model(solver=getsolver())
        @variable(m, -net.maxflow[i] <= q[i=1:narcs] <= net.maxflow[i])

        # take qp = |q|
        @variable(m, qp)
        @constraint(m, qp .>=  q)
        @constraint(m, qp .>= -q)

        @objective(m, :Min, sum(k1*qp + k2*q.^2) +
                            dot(mul, net.A*q + f) + τ*dot(net.A*q+f, net.A*q+f))

        # solve problem
        status = solve(m)

        qa[t, :, ni] = getvalue(q)

        # update multiplier
        pcost += getobjectivevalue(m)
    end

    net.dualcost = pcost
    # update q
    net.q = qa
end


"""Solve transport problem in primal."""
function qsolve!(net::Network)
    k1, k2 = MPTS.Configuration.κ1, MPTS.Configuration.κ2
    k2 *= k1
    narcs = net.narcs
    qa = zeros(Float64, net.ntime-1, size(net.A, 2), net.ninfo)
    pcost = 0.::Float64
    B = net.A[1:end-1, :]'

    for t in 1:net.ntime-1, ni in 1:net.ninfo
        f = @view net.flow[:, t, ni]
        # TODO: dry definition of problem to use hotstart
        m = Model(solver=getsolver())
        @variable(m, -30net.maxflow[i] <= q[i=1:narcs] <= 30net.maxflow[i])

        # equilibrium
        eq = @constraint(m, net.A*q + f .== 0)

        # take qp = |q|
        @variable(m, qp)
        @constraint(m, qp .>=  q)
        @constraint(m, qp .>= -q)

        @objective(m, :Min, sum(k1*qp + k2*q.^2))

        # solve problem
        status = solve(m)
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
    net.primalcost = pcost
    net.q = qa
end


function proj(x, net)

    #= pinv(A)*x =#
    #= A*inv(A'*A)*A'*x =#
    m = JuMP.Model(solver=getsolver())
    @variable(m, -net.maxflow[i] <= q[i=1:narcs] <= net.maxflow[i])

end



"""Solve transport problem in dual."""
function prsolve!(net::Network)

    # we compute analytically the optimal flow
    k1, k2 = MPTS.Configuration.κ1, MPTS.Configuration.κ2
    k2 *= k1

    tcost = 0.
    #= λ = reshape(net.λ, net.ntime-1, size(net.A, 1), net.ninfo) =#
    λ = net.λ
    flows = zeros(Float64, net.ntime-1, size(net.A, 2), net.ninfo)

    for t in 1:net.ntime-1, ni in 1:net.ninfo

        # TODO: dry definition of problem to use hotstart
        m = Model(solver=GurobiSolver(OutputFlag=0))
        mul = @view λ[:, t, ni]
        #= @variable(m, q[i=1:net.narcs]) =#
        @variable(m, -net.maxflow[i] <= q[i=1:net.narcs] <= net.maxflow[i])

        # take qp = |q|
        @variable(m, qp)
        @constraint(m, qp .>=  q)
        @constraint(m, qp .>= -q)

        @objective(m, :Min, sum(k1*qp + k2*q.^2) + dot(mul, net.A*q))

        # solve problem
        solve(m)
        qa = getvalue(q)
        tcost += getobjectivevalue(m)
        flows[t, :, ni] = qa
    end

    #TODO: mem alloc
    net.q = flows
    net.dualcost = -tcost
end

