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
