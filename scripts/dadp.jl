
push!(LOAD_PATH, "..")


using District, StochDynamicProgramming
using Lbfgsb

include("utils.jl")

struct Grid
    ntime
    nodes::Vector{District.AbstractNode}
end


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
end
function DADP(pb::Grid; nsimu=100, nit=10, algo=SDDP(20))

    nnodes = length(pb.nodes)
    ntime = District.ntimesteps(pb.nodes[1].time)


    λ = zeros(Float64, nnodes, ntime-1)
    F = zeros(Float64, nnodes, ntime-1)
    Q = zeros(0)
    scen = [District.genscen(d.model, nsimu) for d in pb.nodes]

    DADP(ntime, Inf, λ, F, Q, algo, nothing, scen, nsimu, nit)
end

function build!(grid::Grid)
    price = zeros(Float64, grid.ntime-1)
    x0 = [.55, 2., 16., 16.]

    conn = PriceInterface(price, GraphConnection(1.))
    for d in pb.nodes
        set!(d, conn)
        District.build!(d, x0)
    end
end



# update graph exchange in nodes subproblems
function swap!(pb::Grid, mul)
    for (id, d) in enumerate(pb.nodes)
        swap!(d, mul[id, :])
    end
end

# solve each node subproblems
function solve!(pb::Grid, dadp::DADP)
    for d in pb.nodes
        solve(d, dadp.algo)
    end
end

function simulate!(pb, dadp)
    dadp.cost = 0
    for d in pb.nodes
        c, f = simulate(d, dadp.scen[id])
        dadp.F[id, :] = f
        dadp.cost = c
    end
end

function setsubgradient!(storage, dadp)
    storage[:] = dadp.F[1, :] - dadp.F[2, :]
end


function oracle(pb::Grid, dadp::DADP)
    f(λ) = dadp.cost

    function grad!(λ, storage)
        # update multiplier
        mul = λ

        swap!(pb, mul)
        solve!(pb, dadp)
        F = simulate!(pb, dadp.nsimu)

        setsubgradient!(storage, pb)
    end

    return f, grad!
end

# Build problem
ts = TimeSpan(200, 3)

# we build two houses
h1 = buildelechouse(ts)
h2 = buildelechouse(ts)

pb = Grid(District.ntimesteps(ts), [h1, h2])
#= algo = DADP(pb) =#


#= f, grad! = oracle(pb, dadp) =#

# Launch Gradient Descent!
    #= return @time Optim.optimize(f, grad!, x0, dadp.algo, options) =#
#= @time lbfgsb(f, grad!, x0; iprint=1, pgtol=1e-8) =#
