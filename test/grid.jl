################################################################################
# District.jl
################################################################################
# Test src/models/grid.jl
################################################################################
push!(LOAD_PATH, "..")

using Base.Test
using District, StochDynamicProgramming, Scenarios

srand(2713)

@testset "Network" begin
    # incidence matrix
    A = [1. -1.]'
    # Time span
    ts = TimeSpan(180, 1)
    nt = District.ntimesteps(ts) -1

    # test wrong definition of network
    @test_throws Exception Network(ts, A, qmax=:Foo)
    @test_throws Exception Network(ts, A, qmax=zeros(10))
    # define network
    net = Network(ts, A, qmax=zeros(1))
    net = Network(ts, A, qmax=6.)
    @test any(x->x==6., net.maxflow)

    nx = District.nnodes(net)
    @test nx == 2

    f = District.flowallocation(net)
    @test any(x->x==0., f)

    @testset "Resolution" begin
        District.swap!(net, rand(nt*nx))
        District.solve!(net)
        @test net.cost ≈ -167.2217 atol=1e-4
        # If ∑ F^i != 0, problem is unfeasible
        @test_throws Exception District.qsolve!(net)
        f = rand(nt)
        ftot = [f; -f]
        District.swap!(net, ftot)
        District.solve!(net)
        @test net.cost ≈ -555.8861 atol=1e-4
    end
end

@testset "Building" begin
    A = [1. -1.]'
    # Time span
    ts = TimeSpan(180, 1)
    nbins = 10
    # we build two houses
    h1 = load(ts, ElecHouse(pv=4, heat=6, bat="bat0", nbins=nbins))
    h2 = load(ts, ElecHouse(pv=0, heat=6, bat="", idhouse=2, nbins=nbins))
    # define initial position
    xini = Dict(h1=> [.55, 2., 20., 20.],
                h2=> [2., 20., 20.])

    # add network
    net = Network(ts, A)

    # Define corresponding grid
    pb = Grid(ts, [h1, h2], net)

    @test District.nnodes(pb) == 2
    @test District.narcs(pb) == 1
    @test District.ntimes(pb) == District.ntimesteps(ts)

    ndev1 = length(h1.devices)
    ndev2 = length(h2.devices)
    # Build SP problems in each node
    for Interface in [FlowInterface, PriceInterface, QuadInterface]
        build!(pb, xini, Interface)
        @test District.checkconsistency(pb, Interface)
        # test that we add only one device corresponding to
        # graph interface
        @test length(h1.devices) == ndev1 + 1
        @test length(h2.devices) == ndev2 + 1
    end

    @test isa(District.initpos(pb), Vector{Float64})

    @testset "Swap" begin
        nt = District.ntimes(pb) -1
        nx = District.nnodes(pb)
        x0 = rand(Float64, nt*nx)

        District.swap!(pb, x0)
        # test if x0 is correctly set in nodes
        @test h1.conn.values == x0[1:nt]
        @test h2.conn.values == x0[(nt+1):end]

        # test if x0 is correctly set in network
        @test net.λ[:, 1] == x0[1:nt]
        @test net.λ[:, 2] == x0[(nt+1):end]
    end

    @testset "Global problem" begin
        # build global sp model
        for nbins in [1, 10]
            sp = District.getproblem(pb, DiscreteLawSampler(nbins, 5, 1))
            @test isa(sp, StochDynamicProgramming.SPModel)
            @test sp.noises[1].supportSize == nbins
        end
    end

    @testset "Simulation" begin
        nassess = 1
        nbins = 20
        sim = Simulator(pb, nassess, sampler=DiscreteLawSampler(nbins, 5, 1))
        @test isa(sim, Simulator)
        @test sim.model.noises[1].supportSize == nbins

        for Samples in [OutSampleScenarios, InSampleScenarios]
            sim = Simulator(pb, nassess, sampler=DiscreteLawSampler(nbins, 5, 1), gen=Samples())
            @test isa(sim, Simulator)
            @test sim.model.noises[1].supportSize == nbins
        end
    end
end
