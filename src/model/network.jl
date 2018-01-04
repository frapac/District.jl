################################################################################
# District.jl
################################################################################
# Electrical networks utilities.
################################################################################
# @TODO/ define properly graph of network (random)


struct Network
    # number of timesteps
    ntime::Int
    # incidence matrix
    A::Array{Float64, 2}
    # number of arcs in graph
    narcs::Int
    # maxflow through arcs
    maxflow::Vector{Float64}
end
function Network(ntime, connex)
    # get number of nodes
    nnodes = size(connex, 1)
    # get number of arcs
    narcs = floor(Int, sum(connex .> 0.)/2)

    A, mf = buildincidence(connex)

    Network(ntime, A, narcs, mf)
end


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
