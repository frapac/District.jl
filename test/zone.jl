################################################################################
# District.jl
################################################################################
# Test src/models/zone.jl
################################################################################

push!(LOAD_PATH, "..")

using Base.Test
using District, StochDynamicProgramming

srand(2713)

nbins = 1
# We build three houses
pb, xini = twelvehouse(nbins=nbins)

nnodes = length(pb.nodes)

# The first two houses are in zone 1, the third in zone 2
membership = vec([1 1 1 3 3 3 3 3 3 2 2 2])
# We fill zone 1 and zone 2
pbreduced = District.reducegrid(pb, membership)
zones = pbreduced.nodes

@testset "Zone" begin
    nbins = 1
    # We build three houses
    pb, xini = twelvehouse(nbins=nbins)

    nnodes = length(pb.nodes)

    # The first two houses are in zone 1, the third in zone 2
    membership = vec([1 1 1 3 3 3 3 3 3 2 2 2])
    # We fill zone 1 and zone 2
    pbreduced = District.reducegrid(pb, membership)
    zones = pbreduced.nodes
    # Test zone assignment
    # Number of zones
    @test length(zones) == 3
    # Number of nodes in zone
    @test length(zones[1].nodes) == 3 && length(zones[2].nodes) == 6 && length(zones[3].nodes) == 3
    # Number of border nodes in zone
    @test length(zones[1].borderindex) == 2 && length(zones[2].borderindex) == 3 && length(zones[3].borderindex) == 1
    # Number of arcs in zones
    @test size(zones[1].net.A,2) == 3 && size(zones[2].net.A,2) == 7 && size(zones[3].net.A,2) == 2
    # Number of arcs in zonal grid
    @test size(pbreduced.net.A, 2) == 4

end

