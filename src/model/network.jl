################################################################################
# District.jl
################################################################################
# Electrical networks utilities.
# - Network implements transportation structure
################################################################################
# @TODO/ define properly graph of network (random)

export Network, Grid, build!

"""
    AbstractNetwork

A Network connects different nodes together.
This object stores the graph's structure.
"""
abstract type AbstractNetwork end

immutable NoneNetwork <: AbstractNetwork end

# TODO: split graph structure from resolution results
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
function Network(ts, A; qmax=6., k1=0., k2=1e-2)
    ntime = ntimesteps(ts)
    nnodes, narcs = size(A)
    λ = zeros(Float64, ntime-1, nnodes)
    F = zeros(Float64, ntime-1, nnodes)
    Q = zeros(Float64, ntime-1, narcs)
    # Build max flow through arcs
    if isa(qmax, Float64)
        maxflow = qmax * ones(Float64, narcs)
    elseif isa(qmax, Vector{Float64}) && (length(qmax) == narcs)
        maxflow = copy(qmax)
    else
        error("Unvalid `qmax`: must be a Float64 or a Vector with size `narcs`")
    end

    return Network(ntime, Inf, Q, A, narcs, maxflow, λ, F, k1, k2)
end

nnodes(net::Network) = size(net.A, 1)

"Set multipliers inside Network `net`."
swap!(net::Network, mul) = net.λ[:] = mul

"""Build incidence matrix and return flows' bounds."""
function buildincidence(connexion::Array{Float64})
    nnodes = size(connexion, 1)
    narcs = floor(Int, sum(connexion .> 0.)/2)

    # incidence matrix
    A = zeros(Int64, nnodes, narcs)

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

"Get laplacian matrix of node-arc incidence matrix `A` with weight q."
function getlaplacian(A::Array{Float64, 2}, q::Array{Float64, 1})
    return A * diagm(q) * A'
end

"Get adjacence matrix of node-arc incidence matrix `A`."
function getadjacence(A::Array{Float64, 2})
    nnodes = size(A, 1)
    L = A * A'
    for i in 1:nnodes
        L[i, i] = 0
    end
    return -L
end

getmaxflow(pos)=sum(build_graph()[pos, pos], 1)[:]

"""
    flowallocation(net::Network)

Compute the injection flows at nodes, as

`` A * q ``

with `A` node-arc incidence matrix of network and `q` current flows
through edges (stored in `net.Q`).
"""
function flowallocation(net::Network)
    dg = zeros(Float64, net.ntime-1, nnodes(net))
    for t in 1:(net.ntime - 1)
        q = net.Q[t, :]
        dg[t, :] = (net.A*q)[:]
    end
    return dg[:]
end

getcost(net::Network, q::Vector{Float64}) = net.k1*sum(abs.(q)) + net.k2*dot(q, q)
