################################################################################
# District.jl
################################################################################
# Test District/simulation
################################################################################

push!(LOAD_PATH, "..")

using Base.Test
using District, StochDynamicProgramming


function buildhouse()

    ts = TimeSpan(0, 1)
    house = House(ts)

    # add devices to house
    for (Device, name) in zip([Battery, ElecHotWaterTank, R6C2], ["bat0", "ehwt0", "rt1988"])
        d = Device(name)
        add!(house, d)
    end

    # add elec load
    wdem = Demands(10, 1)
    add!(house, wdem)

    set!(house, EDFPrice(ts))
    set!(house, ComfortPrice(ts))

    dynam = District.builddynamic(house)
    load = District.buildload(house)

    # add initial pos
    x0 = [.55, 2., 16., 16.]
    District.build!(house, x0)

    return house
end


@testset "Simulation" begin
    house = buildhouse()
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

    for pol in [mpcpol, dppol]
        res = District.simulate(sim, pol)
        @test isa(res, SimulationResult)
    end
end
