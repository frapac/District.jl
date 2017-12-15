
push!(LOAD_PATH, "..")

using Base.Test

using District


@testset "Devices" begin
    for (Stock, ids) in zip([Battery, HotWaterTank, CHP, R6C2],
                            ["bat0", "ehwt0", "chp0", "rt1988"])
        st = Stock(ids)
        @test isa(st, Stock)
    end
end
