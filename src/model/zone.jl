################################################################################
# District.jl
################################################################################
# Generic models for  zone.
################################################################################


mutable struct Zone <: AbstractNodalGrid
    # name of Zone
    name::Symbol
    # Time span
    ts::AbstractTimeSpan
    # Nodes
    nodes::Vector{AbstractNode}
    # Border nodes index
    borderindex::Vector{Int64}
    # Edges
    net::AbstractNetwork
    # Connection interface
    conn::AbstractInterface
    # Bounds on import of the zone
    bounds::Vector{Tuple{Int64,Int64}}
    # SP model
    model
end
function Zone(ts::AbstractTimeSpan,
              nodes::Vector{AbstractNode},
              borderindex::Vector{Int64},
              net::AbstractNetwork)
    Zone(gensym(), ts, nodes, borderindex, net, NoneInterface(), Tuple{Int64,Int64}[], nothing)
end
nbordernodes(pb::Zone) = length(pb.borderindex)
connectionsize(pb::Zone) = (ntimes(pb) - 1) * nbordernodes(pb)
ninjection(pb::Zone) = nbordernodes(pb)
swap!(zone::Zone, exch::Vector{Float64}) = swap!(zone.conn, exch)


###
# Definition of zonal grid

struct ZonalGrid <: AbstractGrid
    # Time span
    ts::AbstractTimeSpan
    # Nodes
    nodes::Vector{Zone}
    # Edges
    net::AbstractNetwork
end
nnodes(pb::ZonalGrid) = sum(nbordernodes.([zone for zone in pb.nodes]))

function build!(grid::ZonalGrid, xini::Dict, Interface::Type=ZoneInterface;
                maxflow=6., tau=1., generation="reduction", nbins=10)
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

"Turns a nodal grid to a zonal grid given a clustering of the nodes"
function reducegrid(pb::Grid, membership::Vector{Int})
    # build zones corresponding to membership
    zones = getzones(pb, membership)
    # build reduced network corresponding to links between zones
    incidence = reducenetwork(pb, zones, membership)
    # update
    updatebounds!(pb, zones, incidence)

    return ZonalGrid(pb.ts, zones, Network(pb.ts, incidence))
end

"Returns the vector containing the zones"
function getzones(pb::Grid, membership::Vector{Int})
    # Zones vector
    zones = Zone[]
    # Adjacency matrix
    adjacencymatrix = getadjacence(pb.net.A)
    # Total number of zones
    nzones = maximum(membership)
    zindex = getzoneorder(membership)

    # Fill the zones vector
    for z in zindex
        # Index of the zone nodes
        belongtozone = membership .== z
        # Extract zone nodes
        zonenode = pb.nodes[belongtozone]
        # Index of border nodes in zone
        borderindex = getborderindex(belongtozone, adjacencymatrix[belongtozone,:])
        # Adjacency matrix of the zone
        zoneadjmat = adjacencymatrix[ belongtozone, belongtozone ]
        # Build incidence matrix of the zone
        A, _ = buildincidence(zoneadjmat)

        push!(zones, Zone(pb.ts, zonenode, borderindex, Network(pb.ts, A)))
    end

    return zones
end

function getzoneorder(membership::Vector{Int})
    zindex = Int64[]
    int=0
    for i in membership
        if i != int
            push!(zindex,i)
            int = i
        end
    end
    return zindex
end

"Returns the vector of indices of the border nodes in a zone"
function getborderindex(belongtozone::BitArray{1}, adjacencymatrix::Array{Float64})
    # number of nodes in zone
    nnodes = sum(belongtozone)
    # Index of border nodes in zone
    borderindex = Int64[]

    for i in 1:nnodes
        # neighbors of i
        neigh = find(x->(x!=0), adjacencymatrix[i, :])
        # checking if a neighbour is out of the zone
        for j in neigh
            if !belongtozone[j]
                push!(borderindex,i)
                break
            end
        end
    end

    return borderindex
end

function reducenetwork(pb::Grid, zones::Vector{Zone}, membership::Vector{Int})
    # Adjacency matrix
    adjacencymatrix = getadjacence(pb.net.A)
    # Index of border nodes
    indexbordernodes = Int64[]
    lastzoneindex=0

    for zone in zones
        indexbordernodes = vcat(indexbordernodes, lastzoneindex + zone.borderindex)
        lastzoneindex += nnodes(zone)
    end
    connexion = adjacencymatrix[ indexbordernodes , indexbordernodes ]
    incidence = buildborderincidence(indexbordernodes, connexion, membership)

    return incidence
end

function buildborderincidence(indexbordernodes::Array{Int64},
                              connexion::Array{Float64}, membership::Vector{Int})
    nbordernodes = length(indexbordernodes)
    # Reduced incidence matrix
    incidence = Float64[]

    # /!\ : We cannot use the function buildincidence because we don't want
    # edges between the border nodes of a same zone to appear
    for i in 1:nbordernodes-1
        for j in i+1:nbordernodes
            if (membership[ indexbordernodes[i] ] != membership[ indexbordernodes[j] ]) &&  (connexion[i, j] == 1)
                vec = zeros(Float64,1,nbordernodes)
                vec[i] = 1
                vec[j] = -1
                incidence = vcat(incidence,vec)
            end
        end
    end
    return Array(incidence')
end

function updatebounds!(pb::Grid, zones::Vector{Zone},incidence::Array{Float64})
    lastzoneindex=0
    for zone in zones
        bounds = Tuple{Int64,Int64}[]

        for i in 1:nbordernodes(zone)
            #number of neighbours out of its zone
            m = sum(abs.(incidence[lastzoneindex + i,:]))
            maxf = m*pb.net.maxflow[1]
            push!(bounds, (-maxf,maxf))
        end

        zone.bounds = bounds
        lastzoneindex+=nbordernodes(zone)
    end
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

            uinj = u[end-ninj+1:end]
            vals = vcat(vals, uinj)
            coefs = vcat(coefs, [ pb.conn.values[t+nb*(ntimes(pb)-1)] for nb in 0:ninj-1 ])
            expr = JuMP.AffExpr(vals, coefs, 0.0)
        end
        return expr
    end

    return costm
end
