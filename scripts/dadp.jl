
push!(LOAD_PATH, "..")


using District, StochDynamicProgramming
using Lbfgsb

include("utils.jl")

struct Grid
    ntime
    nodes::Vector{District.AbstractNode}
end

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

function DADP(pb::Grid; nsimu=100, nit=10, algo=SDDP(10))
    nnodes = length(pb.nodes)
    ntime = District.ntimesteps(pb.nodes[1].time)

    λ = zeros(Float64, nnodes, ntime-1)
    F = zeros(Float64, nnodes, ntime-1)
    Q = zeros(0)
    scen = [District.genscen(d.model, nsimu) for d in pb.nodes]
    mod = Dict()

    DADP(ntime, Inf, λ, F, Q, algo, nothing, scen, nsimu, nit, mod)
end

function build!(grid::Grid)
    price = zeros(Float64, grid.ntime-1)
    x0 = [.55, 2., 16., 16.]

    for d in pb.nodes
        # need one instance of conn per nodes
        # TODO: build! after conn is added
        conn = PriceInterface(copy(price), GraphConnection(1.))
        set!(d, conn)
        District.build!(d, x0)
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
end

# solve each node subproblems
function solve!(pb::Grid, dadp::DADP)
    for d in pb.nodes
        dadp.models[d.name] = solve(d, dadp.algo)
    end
end

function simulate!(pb, dadp)
    dadp.cost = 0
    for (id, d) in enumerate(pb.nodes)
        c, _, u = simulate(dadp.models[d.name], dadp.scen[id])
        dadp.F[id, :] = mean(u[:, :, end], 2)
        dadp.cost += mean(c)
    end
end

function setsubgradient!(storage, dadp)
    storage[:] = dadp.F[1, :] + dadp.F[2, :]
end

∇g(dadp) = (A*(dadp.F[1, :] + dadp.F[2, :])')'[:]


function oracle(pb::Grid, dadp::DADP)
    f(λ) = dadp.cost

    function grad!(λ, storage)
        # update multiplier
        swap!(pb, λ)
        solve!(pb, dadp)
        simulate!(pb, dadp)

        copy!(storage, ∇g(dadp))
    end

    return f, grad!
end

# Build problem
ts = TimeSpan(200, 1)

# we build two houses
h1 = buildelechouse(ts)
h2 = buildelechouse(ts)

pb = Grid(District.ntimesteps(ts), [h1, h2])
build!(pb)
algo = DADP(pb)
solve!(pb, algo)


f, grad! = oracle(pb, algo)

# Launch Gradient Descent!
    #= return @time Optim.optimize(f, grad!, x0, dadp.algo, options) =#
x0 = zeros(Float64, 190)
@time lbfgsb(f, grad!, x0; iprint=1, pgtol=1e-8)
