

struct Grid
    ntime::Int
    # global problem
    model::SPModel
    # production problem (= node)
    nodes::Vector{Node}
    # transport problem (= arcs)
    network::Network
    # number of nodes in problem
    nnodes::Int
    # number of information values
    ninfo::Int
    # DP solver
    Δx
    # import's bound
    fexch
    # decomposition scheme
    dec::Decomposition
end
function Grid(ntime::Int, names::Vector{Symbol}, δx;
              shape=:quad,  # shape of cost's function
              nbins::Int=5, # size of discrete probability distribution
              decom::Decomposition=PriceDecomposition(), # decomposition scheme
             )
    # check consistency of arguments
    @assert length(δx) == length(names)

    # FIXME: generic definition of ny
    ny = 1
    fexch = getmaxflow(getpos.(names))
    # adapt max flow for inner nodes
    fexchn = (isa(decom, PriceDecomposition))?fexch:[Inf for n in names]
    nodes = Node[Node(name, [δx[i]], fexch=fexchn[i], shape=shape, ntime=ntime, nbins=nbins) for (i, name) in enumerate(names)]
    net = Network(ntime, names, ny)

    sp = buildglobal(names, shape=shape, ntime=ntime, nbins=nbins)

    Grid(ntime, sp, nodes, net, length(names), ny, δx, fexch, decom)
end


setexch!(pb::Grid, exch) = setexch!(pb.dec, pb, exch)
setexch!(::PriceDecomposition, pb::Grid, λ)= setmultiplier!(pb, λ)
setexch!(::QuantDecomposition, pb::Grid, F)= setflows!(pb, F)

function setmultiplier!(pb::Grid, λ)
    # update multipliers for production problem
    for (i, node) in enumerate(pb.nodes)
        # must return a view
        #= node.λ'[:] = getmultnode(λ, i, (node.ntime-1)*pb.ninfo) =#
        node.λ = λ[i, :, :]
    end

    # update multipliers for transport problem
    pb.network.λ = λ
end


function setflows!(pb::Grid, flows)
    # update multipliers for production problem
    for (i, node) in enumerate(pb.nodes)
        node.flow = flows[i, :, :]
    end

    # update multipliers for transport problem
    pb.network.flow = flows
end

function getflows(pb::Grid)
    nnodes, ntime, ninfo =  getsize(pb)
    F = zeros(Float64, nnodes, ntime-1, ninfo)
    for (i, node) in enumerate(pb.nodes)
        F[i, :, :,] = node.flow
    end
    F
end

### Dual costs
"""Return dual cost of `Grid` problem."""
dualcost(pb::Grid)=sum(dcost.(pb.nodes)) + dcost(pb.network)

dcost(net::Network) = net.dualcost
dcost(node::Node)   = node.dualcost

### Primal costs
primcost(pb::Grid)=sum(pcost.(pb.nodes)) + pcost(pb.network)

pcost(net::Network) = net.primalcost
pcost(node::Node)   = node.primalcost

# PROJECTION
function projection!(x, pb::Grid)
    net = pb.network

    for t in 1:net.ntime-1, ni in 1:net.ninfo
        f = x[:, t, ni]
        # projection onto ∑ x_i = 0 subspace
        x[:, t, ni] = proj(f)
    end
end
function proj(x)
    a = ones(length(x))
    λ = a'*x/(a'*a)
    x - λ*a
end
function proj2(x)
    imin = indmin(x)
    rest = 0.
    for i in 1:endof(x)
        if i != imin
            rest += x[i]
        end
    end
    x[imin] = -rest
    x
end

# UTILS
getnames(pt::Grid)=Symbol[node.name for node in pt.nodes]

getsize(pb::Grid) = (pb.nnodes, pb.ntime, pb.ninfo)
