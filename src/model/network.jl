################################################################################
# District.jl
################################################################################
# Electrical networks utilities.
# - Network implements transportation structure
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
    # node imports
    F::Array{Float64, 2}
    # transport cost
    k1::Float64
    k2::Float64
end
function Network(ts, A)
    ntime = ntimesteps(ts)
    nnodes, narcs = size(A)
    λ = zeros(Float64, ntime-1, nnodes)
    F = zeros(Float64, ntime-1, nnodes)
    Q = zeros(Float64, ntime-1, narcs)
    # TODO: dry
    maxflow = 6 * ones(Float64, narcs)

    return Network(ntime, Inf, Q, A, narcs, maxflow, λ, F, 0., 1e-2)
end

nnodes(net::Network) = size(net.A, 1)

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

function flowallocation(net::Network)
    dg = zeros(Float64, net.ntime-1, nnodes(net))
    for t in 1:(net.ntime - 1)
        q = net.Q[t, :]
        dg[t, :] = (net.A*q)[:]
    end
    return dg[:]
end
