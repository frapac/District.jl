################################################################################
# District.jl
################################################################################
# Generic models for  zone.
################################################################################

include("../../scripts/graph.jl")

mutable struct Zone <: AbstractNodalGrid
    # name of Zone
    name::Symbol
    # Time span
    ts::AbstractTimeSpan
    # Nodes
    nodes::Vector{AbstractNode}
    # Border nodes index
    borderindex::Dict{AbstractNode, Int64}
    # Edges
    net::AbstractNetwork
    # Connection interface
    conn::AbstractInterface
    # SP model
    model
end

struct ZonalGrid <: AbstractGrid
    # Time span
    ts::AbstractTimeSpan
    # Nodes
    nodes::Vector{Zone}
    # Edges
    net::AbstractNetwork
end

function Zone(ts::AbstractTimeSpan, 
    nodes::Vector{AbstractNode}, 
    borderindex::Dict{AbstractNode, Int64}, 
    net::AbstractNetwork)
    Zone(gensym(), ts, nodes, borderindex, net, NoneInterface(),nothing)
end

nnodes(pb::ZonalGrid) = sum(length.([zone.borderindex for zone in pb.nodes]))
nbordernodes(pb::Zone) = length(pb.borderindex)
connectionsize(pb::Zone) = (ntimes(pb) - 1) * nbordernodes(pb)
ninjection(pb::Zone) = nbordernodes(pb)


function swap!(zone::Zone, exch::Vector{Float64})
    swap!(zone.conn, exch)
end

function build!(grid::ZonalGrid, xini::Dict, Interface::Type=ZoneInterface;
                maxflow=6., tau=1., generation="reduction",nbins=10)
    for zone in grid.nodes
        price = zeros(Float64, connectionsize(zone))
        linker = zeros(Int64, nnodes(zone), nbordernodes(zone))

        j = 1
        for i in collect(values(zone.borderindex))
            linker[i, j] = 1
            j += 1
        end

        conn = Interface(price, linker)
        build!(zone, xini, PriceInterface)
        zone.conn = conn
        zone.model = getproblem(zone, generation, nbins)
    end

end


function reducepb(pb::Grid, membership::Vector{Int})
    

    
    # Instantiate zones 
    zones = Zone[]

    fillzones!(pb, zones, membership)

    incidence = reducenetwork(pb, zones, membership)


    return ZonalGrid(pb.ts, zones, Network(pb.ts, incidence)) 
end

function fillzones!(pb::Grid, zones::Vector{Zone}, membership::Vector{Int})
    
    # Adjacency matrix
    adjacencymatrix = getadjacence(pb.net.A)
    # Total number of nodes
    nnodes = size(pb.nodes,1)
    # Total number of zones
    nzones = maximum(membership)
    # Map between nodes and their index
    index = Dict(pb.nodes[i]=> i for i in 1:nnodes)

    # Fill the vectors zonenode, and borderindex
    for z in 1:nzones
        # Vector that associates a zone to its nodes
        zonenode = AbstractNode[]
        # Vector that associates a zone to a mapping between border nodes and their index in the node list 
        borderindex = Dict{AbstractNode,Int64}()
        # Index of a border node in the node list
        bindex = 0

        # fill nodes in zone z
        for i in 1:nnodes
            # If node i belongs to zone z
            if membership[i] == z
                push!(zonenode, pb.nodes[i])
                bindex += 1    
                # If node i has a neighbour out of its zone
                for j in 1:nnodes
                    if (membership[i] != membership[j]) && (adjacencymatrix[i,j] == 1)
                        get!(borderindex, pb.nodes[i], bindex)
                        break
                    end
                end
            end
        end
    

        nzonenodes = size(zonenode,1)
        adjmat = adjacencymatrix[ [index[zonenode[i]] for i in 1:nzonenodes] , [index[zonenode[i]] for i in 1:nzonenodes] ]
        
        A, bounds = buildincidence(adjmat)

        push!(zones, Zone(pb.ts, zonenode, borderindex, Network(pb.ts, A)))
    end
end

function reducenetwork(pb::Grid, zones::Vector{Zone}, membership::Vector{Int})
    
    # Adjacency matrix
    adjacencymatrix = getadjacence(pb.net.A)
    # Total number of nodes
    nnodes = size(pb.nodes,1)
    # Map between nodes and their index
    index = Dict(pb.nodes[i]=> i for i in 1:nnodes)

    # Total number of border nodes
    nbordernodes = sum(length.([zone.borderindex for zone in zones]))
    # Vector of all border nodes
    allbordernodes = AbstractNode[]
    for zone in zones
        for i in collect(values(zone.borderindex))
            push!(allbordernodes, zone.nodes[i])
        end
    end

    connexion = adjacencymatrix[ [index[allbordernodes[i]] for i in 1:nbordernodes] , [index[allbordernodes[i]] for i in 1:nbordernodes] ]
    nborderarcs = floor(Int, sum(connexion .> 0.)/2)
    
    # Reduced incidence matrix 
    incidence = zeros(Float64, nbordernodes, nborderarcs)
    iarc = 1

    for i in 1:nbordernodes-1
        for j in i+1:nbordernodes
            if (membership[index[allbordernodes[i]]] != membership[index[allbordernodes[j]]]) &&  (adjacencymatrix[ index[allbordernodes[i]], index[allbordernodes[j]] ] == 1)
                incidence[i, iarc] = 1
                incidence[j, iarc] = -1
                iarc += 1
            end
        end
    end

    return incidence
end

function objective(pb::Zone)
    ninj = ninjection(pb)

    function costm(m, t, x, u, w)
        vals = JuMP.Variable[]
        coefs = Float64[]

        # add decomposition price
        if isa(pb.conn, ZoneInterface)
            # add < λ, F_out >
            u = m[:u]

            vcat(vals, u[end-ninj+1:end])
            vcat(coefs, [ pb.conn.values[t+nb*(ntimes(pb)-1)] for nb in 0:ninj-1 ])
            expr = JuMP.AffExpr(vals, coefs, 0.0)
            for i in 1:ninj
                expr += JuMP.QuadExpr([u[end-i+1]], [u[end-i+1]], [1e-2], 0.)
            end
        end
        return expr
    end

    return costm
end