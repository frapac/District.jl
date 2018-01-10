################################################################################
# District.jl
################################################################################
# Define whole microgrid with different interconnected noces.
################################################################################

# TODO: interface properly with LightGraphs
# Grid is the metastructure storing all interconnected elements in the grid
struct Grid
    name::Symbol
    # time period considered
    time::AbstractTimeSpan
    # global problem
    model::SPModel
    # production problem (= node)
    nodes::Vector{AbstractNode}
    # transport problem (= arcs)
    network::Network
    # import's bound
    fexch
end
function Grid(ts::AbstractTimeSpan)

end
