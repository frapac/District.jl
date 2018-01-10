################################################################################
# District.jl
################################################################################
# Test District/models/buildings.jl
################################################################################

push!(LOAD_PATH, "..")

using Base.Test

using District


@testset "Building" begin
    ts = TimeSpan(0, 1)
    house = House(ts)
    @test isa(house, House)

    # add three devices to house
    b = Battery("bat0")
    add!(house, b)
    @test length(house.devices) == 1

    hwt = ElecHotWaterTank("ehwt0")
    add!(house, hwt)
    @test nstocks(house) == 2

    thm = R6C2("rt1988")
    add!(house, thm)
    @test nstocks(house) == 4

    wdem = Demands(10, 1)
    add!(house, wdem)
    @test nnoises(house) == 2

    # add elec and comfort price:
    set!(house, EDFPrice(ts))
    set!(house, ComfortPrice(ts))

    dynam = District.builddynamic(house)
    @test isa(dynam(1, [4, 0, 0, 0], zeros(10), zeros(3)), Vector{Float64})

    xb = District.xbounds(house)
    @test length(xb) == nstocks(house)
    District.ubounds(house)

    bcost = District.objective(house)
    @test isa(bcost, Function)

    # add initial pos
    x0 = [.55, 6., 16., 16.]
    District.build!(house, x0)
    @test isa(house.model, District.LinearSPModel)
end
