
push!(LOAD_PATH, "..")

using Base.Test

using District


@testset "Devices" begin

    bat = Battery("bat0")
    @test isa(bat, Battery)

    hwt = HotWaterTank("ehwt0")
    @test isa(hwt, HotWaterTank)

    chp = CHP("chp0")
    @test isa(chp, CHP)


    th  = R6C2("rt1988")
end

