
push!(LOAD_PATH, "..")

using Base.Test

using District


@testset "Building" begin
    ts = TimeSpan(0, 1)
    house = House(ts)
    @test isa(house, House)

    b = Battery("bat0")
    add!(house, b)
    @test length(house.devices) == 1
    ex = District.parsedevice(b, 1, 1, .25)

    hwt = HotWaterTank("ehwt0")
    add!(house, hwt)
    ex = District.parsedevice(hwt, 1, 1, .25)

    thm = R6C2("rt1988")
    add!(house, thm)

    wdem = Demands()
    add!(house, wdem)

    dynam = District.builddynamic(house)
    @test isa(dynam(1, [4, 0, 0, 0], zeros(10), zeros(3)), Vector{Float64})
    load = District.buildload(house)
    @test isa(load(40, [4, 0, 0, 0], zeros(10), zeros(3)), Float64)

    xb = District.xbounds(house)
    @test length(xb) == nstocks(house)
    District.ubounds(house)

    bcost = District.objective(house)

    x0 = [.55, 6., 16., 16.]
    District.build!(house, 10, x0)
end
