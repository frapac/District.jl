################################################################################
# District.jl
################################################################################
# Generic models for  zone.
################################################################################

include("../../scripts/graph.jl")

mutable struct Zone <: AbstractNodalGrid
    # name of House
    name::Symbol
    # Time span
    ts::AbstractTimeSpan
    # Nodes
    nodes::Vector{AbstractNode}
    # Nodes
    bordernodes::Vector{AbstractNode}
    # Border nodes index
    borderindex::Dict{AbstractNode, Int64}
    # Edges
    net::AbstractNetwork
    # Inteface with graph
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

Zone(ts::AbstractTimeSpan, nodes::Vector{AbstractNode}, bordernodes::Vector{AbstractNode}, borderindex::Dict{AbstractNode, Int64}, net::AbstractNetwork)= Zone(gensym(), ts, 
                                                                                                                                                        nodes, bordernodes, 
                                                                                                                                                        borderindex, net, 
                                                                                                                                                        NoneInterface(),nothing)

nnodes(pb::ZonalGrid) = sum(length.([zone.bordernodes for zone in pb.nodes]))
connectionsize(zone::Zone) = (ntimes(zone) - 1) * length(zone.bordernodes)

function swap!(zone::Zone, exch::Vector{Float64})
    slambda = 0
    for d in zone.bordernodes
        s1 = slambda + 1
        slambda += connectionsize(d)
        s2 = slambda
        swap!(d, exch[s1:s2])
    end
end

function zonebuild!(grid::ZonalGrid, xini::Dict, Interface::Type=PriceInterface;
                maxflow=6., tau=1.)
    for zone in grid.nodes
        build!(zone, xini, PriceInterface)
        zone.model = getproblem(zone, "reduction", 30)
    end

end


function build!(zone::Zone, xini::Dict, Interface::Type=PriceInterface;
                maxflow=6., tau=1.)
    conn = NoneInterface()
    
    for d in zone.nodes
        price = zeros(Float64, ntimes(zone)-1)
        # instantiate connection interface
        conn = Interface(price, GraphConnection(maxflow))
        
        # add connection to particular device
        set!(d, conn)

        build!(d, xini[d])
    end

    for d in zone.bordernodes
        price = zeros(Float64, ntimes(zone)-1)
        # instantiate connection interface
        conn = Interface(price, InOutGraphConnection(maxflow))
        
        # add connection to particular device
        set!(d, conn)

        build!(d, xini[d])
    end

    zone.conn = conn
end



function initvector(vector::Vector{Vector{Float64}})
    for i in 1:size(vector,1)
        vector[i] = []
    end
end

function fillzones(pb::Grid, membership::Vector{Int})
    
    # Total number of zones
    nzones = maximum(membership)
    # Total number of nodes
    nnodes = size(pb.nodes,1)

    # Adjacency matrix
    adjmatrix = getadjacence(pb.net.A)
    
    # Instantiate zones 
    zones = Zone[]

    # Map between nodes and their index
    index = Dict(pb.nodes[i]=> i for i in 1:nnodes)

    # Fill the vectors zonenode, zoneborder and borderindex
    for z in 1:nzones
         # Vector that associates a zone to its nodes
        zonenode = AbstractNode[]
        # Vector that associates a zone to its border nodes
        zoneborder = AbstractNode[]
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
                    if (membership[i] != membership[j]) && (adjmatrix[i,j] == 1)
                        push!(zoneborder, pb.nodes[i])
                        get!(borderindex, pb.nodes[i], bindex)
                        break
                    end
                end
            end
        end
    

        nznodes = size(zonenode,1)
        nzarcs = sum(sum(adjmatrix[ [index[zonenode[i]] for i in 1:nznodes] , [index[zonenode[i]] for i in 1:nznodes] ]))/2
        incidence = zeros(Float64, nznodes, nzarcs)
        iarc = 1

        for i in 1:nznodes-1
            for j in (i+1):nznodes
                if adjmatrix[ index[zonenode[i]] , index[zonenode[j]] ]  == 1 
                    incidence[i, iarc] = 1
                    incidence[j, iarc] = -1
                    iarc += 1
                end
            end
        end

        push!(zones, Zone(pb.ts, zonenode, zoneborder, borderindex, Network(pb.ts, incidence)))
    end

    # Total number of border nodes
    nbordernodes = sum(length.([zone.bordernodes for zone in zones]))
    # Vector of all border nodes
    allbordernodes = AbstractNode[]
    for zone in zones
        for h in zone.bordernodes
            push!(allbordernodes, h)
        end
    end
    nborderarcs = sum(sum(adjmatrix[ [index[allbordernodes[i]] for i in 1:nbordernodes] , [index[allbordernodes[i]] for i in 1:nbordernodes] ]))/2

    # Reduced incidence matrix 
    incidence = zeros(Float64, nbordernodes, nborderarcs)
    iarc = 1

    for i in 1:nbordernodes-1
        for j in i+1:nbordernodes
            if (membership[index[allbordernodes[i]]] != membership[index[allbordernodes[j]]]) &&  (adjmatrix[ index[allbordernodes[i]], index[allbordernodes[j]] ] == 1)
                incidence[i, iarc] = 1
                incidence[j, iarc] = -1
                iarc += 1
            end
        end
    end


    return zones, Network(pb.ts, incidence)
end
