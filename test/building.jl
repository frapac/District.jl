
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
end
