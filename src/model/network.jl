################################################################################
# District.jl
################################################################################
# Electrical networks utilities.
# - Network implements transportation structure
# - Grid implements global network (nodes + interconnected edges)
################################################################################
# @TODO/ define properly graph of network (random)

export Network, Grid, build!

abstract type AbstractNetwork end

immutable NoneNetwork <: AbstractNetwork end

mutable struct Network <: AbstractNetwork
    # number of timesteps
    ntime::Int
    # current transport cost
    cost::Float64
    # flow through network
    Q::Array{Float64, 2}
    # incidence matrix
    A::Array{Float64, 2}
    # n arcs
    narcs::Int
    # maximum flow through arcs
    maxflow::Vector{Float64}
    # multiplier
    λ::Array{Float64, 2}
end
function Network(ts, A)
    ntime = ntimesteps(ts)
    nnodes, narcs = size(A)
    λ = zeros(Float64, ntime-1, nnodes)
    Q = zeros(Float64, ntime-1, narcs)
    # TODO: dry
    maxflow = 6 * ones(Float64, narcs)

    return Network(ntime, Inf, Q, A, narcs, maxflow, λ)
end

"Set multipliers inside Network `net`."
swap!(net::Network, mul) = net.λ[:] = mul

"""Build incidence matrix and return flows' bounds."""
function buildincidence{T}(connexion::Matrix{T})
    nnodes = size(connexion, 1)
    narcs = floor(Int, sum(connexion .> 0.)/2)

    # incidence matrix
    A = zeros(Int, nnodes, narcs)

    ic = 0
    bounds = Float64[]
    for ix in 1:(nnodes-1), iy in (ix+1):nnodes
        if connexion[ix, iy] > 0
            ic += 1
            A[ix, ic] =  1
            A[iy, ic] = -1
            push!(bounds, connexion[ix, iy])
        end
    end
    return A, bounds
end

getmaxflow(pos)=sum(build_graph()[pos, pos], 1)[:]



################################################################################
# Define whole microgrid with different interconnected noces.
abstract type AbstractGrid end

struct Grid <: AbstractGrid
    # Time span
    ts::AbstractTimeSpan
    # Nodes
    nodes::Vector{AbstractNode}
    # Edges
    net::AbstractNetwork
end
Grid(ts::AbstractTimeSpan) = Grid(ts, AbstractNode[], NoneNetwork())

nnodes(pb::Grid) = length(pb.nodes)
narcs(pb::Grid) = pb.net.narcs
ntimes(pb::Grid) = ntimesteps(pb.ts)


# Build models inside `grid`
function build!(grid::Grid, xini)
    for d in grid.nodes
        price = zeros(Float64, ntimes(grid)-1)
        # instantiate connection interface
        # TODO: dry connection size
        conn = PriceInterface(price, GraphConnection(6.))
        # add connection to particular device
        set!(d, conn)
        # build SP model inside `d`
        build!(d, xini[d])
    end
end

# update graph exchange in nodes subproblems
function swap!(pb::Grid, mul::Vector{Float64})
    ntime = ntimes(pb) - 1
    for (id, d) in enumerate(pb.nodes)
        s1 = (id-1) * ntime + 1
        s2 = id * ntime
        swap!(d, mul[s1:s2])
    end
    # swap transport problem
    swap!(pb.net, mul)
end

function getproblem(pb::Grid, x0::Vector{Float64})
    xb = vcat(xbounds.(pb.nodes)...)
    ub = vcat(ubounds.(pb.nodes)...)
    # add flow control
    ub = vcat(ub, [(-m, m) for m in pb.net.maxflow])

    dynam = builddynamic(pb)
    costm = buildcost(pb)
    constr = buildconstr(pb)
    laws = buildlaws(pb)

    spmodel = StochDynamicProgramming.LinearSPModel(ntimes(pb), ub,
                                                  x0, costm,
                                                  dynam,
                                                  tonoiselaws(laws),
                                                  info=:HD,
                                                  eqconstr=constr,
                                                  Vfinal=fcost)

    set_state_bounds(spmodel, xb)
    return spmodel
end

function builddynamic(pb::Grid)
    ntime = ntimes(pb)
    xindex = 1
    uindex = 1
    exdyn = Expr(:vect)

    params = Dict()
    params["text"] = loadweather(OutdoorTemperature(), ntime)

    for dev in pb.nodes
        # update irradiation in params
        pint, pext = get_irradiation(dev)
        params["pint"] = pint
        params["pext"] = pext

        # parse dynamics of nodes
        # TODO: w index in builder
        ex = parsebuilding(dev, xindex, uindex, dev.time.δt, params)
        xindex += nstates(dev)
        uindex += ncontrols(dev)
        # push expression in global expression
        push!(exdyn.args, ex...)
    end

    return eval(:((t, x, u, w) -> $exdyn))
end

function buildcost(pb::Grid)
    function costm(m, t, x, u, w)
        cost = AffExpr(0.)
        uindex = 1
        for d in pb.nodes
            # rebuild problem in Node `d`, starting from uindex
            buildload!(d, uindex)
            cost += objective(d)(t, x, u, w)
            uindex += ncontrols(d)
        end
        return cost
    end
end

# get index of flow importations in node
# remember: flow is last control in each node (u[end]) so its position
# is ncontrols
getflowindex(pb::Grid) = cumsum(ncontrols.(pb.nodes))

# Build coupling constraint Aq + F for a grid.
function buildconstr(pb::Grid)
    na = narcs(pb)
    A = pb.net.A
    # TODO: build
    findex = getflowindex(pb)

    constr(t, x, u, w) = A * u[end-na+1:end] + u[findex]
end
