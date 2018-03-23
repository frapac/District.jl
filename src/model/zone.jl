################################################################################
# District.jl
################################################################################
# Generic models for  zone.
################################################################################


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
    borderindex::Dict{AbstractNode, Float64}
    # Edges
    net::AbstractNetwork
    # SP model
    model::StochDynamicProgramming.SPModel
end

struct ZonalGrid <: AbstractGrid
    # Time span
    ts::AbstractTimeSpan
    # Nodes
    nodes::Vector{Zone}
    # Edges
    net::AbstractNetwork
end

Zone(ts::AbstractTimeSpan) = Zone(gensym(), ts, AbstractNode[], AbstractNode[], Dict(), NoneNetwork(), nothing)

nnodes(pb::ZonalGrid) = sum(length.(zone.bordernodes for zone in pb.nodes))
sizelambda(pb::Zone) = ntimes(pb) * size(bordernodes)

function swap!(zone::Zone, exch::Array{Float64})
	for (id,d) in zone.bordernodes
		slambda = sizelambda(d)
        s1 = (id-1) * slambda + 1
        s2 = id * slambda
		swap!(d, exch[s1:s2])
	end
end

function zonebuild!(grid::ZonalGrid, xini::Dict, Interface::Type=PriceInterface;
                maxflow=6., tau=1.)
	for zone in grid.nodes
		build!(zone, xini[d in zone.nodes], PriceInterface)
		zone.model = getproblem(zone, generation="reduction", nbins=30)
	end

end


function build!(zone::Zone, xini::Dict, Interface::Type=PriceInterface;
                maxflow=6., tau=1.)
	

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
end



function initvector(vector::Vector{Vector{Any}})
	for i in 1:size(vector)
		vector[i] = []
	end
end

function fillzones(pb::Grid, membership::Vector{Int})
	
	# Total number of zones
	nzones = maximum(membership)
	# Total number of nodes
	nbnodes = size(pb.nodes)

 	# Adjacency matrix
	adjmatrix = getadjacence(pb.net.A)
	
	# Vector that associates a zone to its nodes
	zonenode = Vector{Vector{AbstractNode}}(nzones)
	# Vector that associates a zone to its border nodes
	zoneborder = Vector{Vector{AbstractNode}}(nzones)
	# Vector that associates a zone to a mapping between border nodes and their index in the node list 
	borderindex = Vector{Dict{AbstractNode,Float64}}(nzones)

	# Initiate vector of vectors with empty vectors
	initvector(zonenode)
	initvector(zoneborder)

	

	
	# Fill the vectors zonenode, zoneborder and borderindex
	for z in 1:nzones
		# Index of a border node in the node list
		bindex = 0.
	    # fill nodes in zone z
	    for i in 1:nbnodes
	        # If node i belongs to zone z
	        if membership[i] == z
	            append!(zonenode[z], pb.nodes[i])
	            bindex += 1    
                # If node i has a neighbour out of its zone
                for j in 1:nbnodes
                    if (membership[i] != membership[j]) && (adjmatrix[i,j] == 1)
                        append!(zoneborder[z], pb.nodes[i])
                        get!(borderindex[z], pb.nodes[i], bindex)
                        break
                    end
                end
	        end
	    end
	end

	# Instantiate zones 
	zones = Vector{Zone}()

	# Map between nodes and their index
	index = Dict(pb.nodes[i]=> i for i in 1:nbnodes)
	
	# Build network of each zone and fill the vector zones
	for z in 1:nzones
		# Build network of zone z
	    znodes = zonenode[z]
	    nznodes = size(znodes)
	    incidence = Vector{Vector{Float64}}(nznodes)
	    initvector(incidence)

	    for i in 1:nznodes-1
			for j in (i+1):nznodes
				if adjmatrix[index[znodes[i]],index[znodes[j]]] == 1
					vect = zeros(Float64,nznodes)
					vect[i] = 1
					vect[j] = -1
					append!(incidence, vect)
				end
			end
		end

	    append!(zones, Zone(ts, zonenode[z], zoneborder[z], borderindex[z], Network(ts, incidence)))
	end


	# Total number of border nodes
	nbordernodes = sum(length.(zone.bordernodes for zone in zones))
	# Vector of all border nodes
	allbordernodes = Vector{AbstractNode}()
	for zone in zones
		for h in zone.bordernodes
			append!(allbordernodes, h)
		end
	end

	# Reduced incidence matrix 
	incidence = Vector{Vector{Float64}}(nbordernodes)
	initvector(incidence)

	for i in 1:nbordernodes-1
		for j in i+1:nbordernodes
			if (membership[index[allbordernodes[i]]] != membership[index[allbordernodes[j]]]) &&  (adjmatrix[ index[allbordernodes[i]], index[allbordernodes[j]] ] == 1)
				vect = zeros(Float64,nbordernodes)
				vect[i] = 1
				vect[j] = -1
				append!(incidence, vect)
			end
		end
	end


	return zones, Network(ts, incidence)
end