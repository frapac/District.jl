
push!(LOAD_PATH, "..")


using District, StochDynamicProgramming
using Lbfgsb, Optim, JuMP, Gurobi

include("network.jl")
srand(2713)

struct Grid
    # TODO: redundant information
    ntime::Int
    ts::District.AbstractTimeSpan
    nodes::Vector{District.AbstractNode}
    net::Network
end

nnodes(pb::Grid) = length(pb.nodes)
narcs(pb::Grid) = pb.net.narcs

A = [1. -1.]'

mutable struct DADP
    ntime::Int
    cost::Float64
    λ::Array{Float64, 2}
    F::Array{Float64, 2}
    Q::Array{Float64}
    algo::District.AbstractSolver
    pb
    scen::Array
    nsimu::Int
    nit::Int
    models::Dict
end

function DADP(pb::Grid; nsimu=100, nit=10, algo=SDDP(nit))
    nnodes = length(pb.nodes)
    ntime = District.ntimesteps(pb.nodes[1].time)

    λ = zeros(Float64, nnodes, ntime-1)
    F = zeros(Float64, nnodes, ntime-1)
    Q = zeros(Float64, ntime-1, pb.net.narcs)
    scen = [District.genscen(d.model, nsimu) for d in pb.nodes]
    mod = Dict()

    DADP(ntime, Inf, λ, F, Q, algo, nothing, scen, nsimu, nit, mod)
end

function build!(grid::Grid, xini)
    price = zeros(Float64, grid.ntime-1)

    for d in pb.nodes
        # need one instance of conn per nodes
        # TODO: build! after conn is added
        conn = PriceInterface(copy(price), GraphConnection(6.))
        set!(d, conn)
        #= set!(d, RecoursePrice(grid.ts)) =#
        District.build!(d, xini[d])
    end
end

# update graph exchange in nodes subproblems
function swap!(pb::Grid, mul)
    ntime = pb.ntime - 1
    for (id, d) in enumerate(pb.nodes)
        s1 = (id-1) * ntime + 1
        s2 = id * ntime
        District.swap!(d, mul[s1:s2])
    end
    # swap transport problem
    swap!(pb.net, mul)
end

# solve each node subproblems
function solve!(pb::Grid, dadp::DADP)
    for d in pb.nodes
        dadp.models[d.name] = District.solve(d, dadp.algo)
    end
    # solve transport problem
    solve!(pb.net)
end

function simulate!(pb, dadp)
    dadp.cost = 0
    for (id, d) in enumerate(pb.nodes)
        c, _, u = simulate(dadp.models[d.name], dadp.scen[id])
        dadp.F[id, :] = mean(u[:, :, end], 2)
        dadp.cost += mean(c)
    end

    dadp.cost += pb.net.cost
    # TODO: not in the same sense
    copy!(dadp.Q, pb.net.Q)
end

function ∇g(pb::Grid, dadp::DADP)
    dg = zeros(Float64, dadp.ntime-1, nnodes(pb))
    for t in 1:(dadp.ntime - 1)
        f = dadp.F[:, t]
        q = dadp.Q[t, :]
        dg[t, :] = (pb.net.A*q + f)[:]
    end
    return dg[:]
end

function oracle(pb::Grid, dadp::DADP)
    # take care: we aim at find a maximum (min f = - max -f )
    f(λ) = - dadp.cost

    function grad!(λ, storage)
        # update multiplier
        swap!(pb, λ)
        solve!(pb, dadp)
        simulate!(pb, dadp)

        copy!(storage, -∇g(pb, dadp))
    end

    return f, grad!
end

# Build problem
ts = TimeSpan(200, 1)

# we build two houses
h1 = load(ts, ElecHouse(pv=4, heat=6, bat="bat0", nbins=1))
h2 = load(ts, ElecHouse(pv=0, heat=6, bat="", idhouse=2, nbins=1))
xini = Dict(h1=> [.55, 2., 20., 20.],
            h2=> [2., 20., 20.])

net = Network(ts, A)
pb = Grid(District.ntimesteps(ts), ts, [h1, h2], net)
build!(pb, xini)
algo = DADP(pb, nsimu=1)
solve!(pb, algo)


f, grad! = oracle(pb, algo)

# Launch Gradient Descent!
x0 = zeros(Float64, 190)
res = @time lbfgsb(f, grad!, x0; iprint=1, pgtol=1e-8)
#= options = Optim.Options(g_tol=1e-8, allow_f_increases=false, show_trace=true, =#
#=                         iterations=20, store_trace=true, extended_trace=false) =#
#= res = @time Optim.optimize(f, grad!, x0, BFGS(), options) =#
