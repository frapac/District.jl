################################################################################
# District.jl
################################################################################
# Generic models for  zone.
################################################################################


mutable struct Zone <: AbstractNode
	# Time span
    ts::AbstractTimeSpan
    # Nodes
    nodes::Vector{AbstractNode}
    # Nodes
    bordernodes::Vector{AbstractNode}
    # Edges
    net::AbstractNetwork
    # SP model
    model
end

Zone(ts::AbstractTimeSpan) = Zone(ts, AbstractNode[], NoneNetwork())

nnodes(pb::Zone) = length(pb.nodes)
narcs(pb::Zone) = pb.net.narcs
ntimes(pb::Zone) = ntimesteps(pb.ts)

sizelambda(pb::Zone) = ntimes(pb) * size(bordernodes)

swap!(zone::AbstractNode, exch::Array{Float64})
	for (id,d) in zone.bordernodes
		slambda = sizelambda(d)
        s1 = (id-1) * slambda + 1
        s2 = id * slambda
		swap!(d, exch[s1:s2])
	end
end

zonebuild!(grid::Grid, xini::Dict, Interface::Type=PriceInterface;
                maxflow=6., tau=1.)
	for zone in grid.nodes
		build!(zone, xini[d in zone.nodes], PriceInterface)
	end
end


build!(zone::AbstractNode, xini::Dict, Interface::Type=PriceInterface;
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