################################################################################
# District.jl
################################################################################
# Test District/simulation
################################################################################

push!(LOAD_PATH, "..")

using Base.Test
using District, StochDynamicProgramming

srand(2713)


@testset "Simulation" begin
    ts = TimeSpan(0, 1)
    house = load(ts, ElecHouse())
    District.build!(house, [2., 16., 16.])
    nscen = 10

    # Build simulator
    sim = Simulator(house, nscen)
    @testset "Simulator" begin
        ntime = District.ntimesteps(house.time)
        @test isa(sim, Simulator)
        @test sim.ts == house.time
        @test size(sim.scenarios) == (ntime, nscen, nnoises(house))
    end

    # Build MPC policies
    forecast = District.genforecast(house.time, house.noises)
    @test isa(forecast, Array{Float64, 2})
    mpcpol = District.MPCPolicy(forecast)

    # Fetch null value functions
    sddp = District.solve(house, SDDP(5))
    V = sddp.bellmanfunctions
    dppol = District.HereAndNowDP(V)

    for (pol, vals) in zip([mpcpol, dppol], [14.2567, 15.4450])
        res = District.simulate(sim, pol)
        @test isa(res, SimulationResult)
        @test .25*mean(res.costs) ≈ vals  atol=1e-4
        println(res)
    end

    # test scenarios generator
    for gen in [InSampleScenarios(), OutSampleScenarios()]
        sim = Simulator(house, nscen, gen)
        @test size(sim.scenarios, 2) == nscen
    end
end
