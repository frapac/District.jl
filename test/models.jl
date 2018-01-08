################################################################################
# District.jl
################################################################################
# Test District/models/
################################################################################

push!(LOAD_PATH, "..")

using Base.Test
using District, Scenarios

@testset "District model" begin
    @testset "Data" begin
        d = District.loadweather()
        @test isa(d, Array{Float64, 2})
        @test size(d) == (35037, 16)
    end


    @testset "Utils" begin
        @testset "TimeSpan" begin
            ts = TimeSpan(1, 7)
            @test District.ntimesteps(ts) == 7 * 96
            @test isa(District.unravel(ts), Tuple{Int64, Int64})
            @test isa(District.weekcycle(ts), Vector{Int64})
            ts = TimeSpan(1, 30)
            @test isa(District.weekcycle(ts), Vector{Int64})
        end
    end


    @testset "Data" begin
        ts = TimeSpan(0, 1)

        @testset "Price" begin
            for tariff in [EDFPrice, EPEXPrice]
                price =  loadprice(tariff(), ts)
                @test isa(price, Array{Float64})
                @test length(price) == 96
            end
        end
        @testset "Setpoint" begin
            stp = loadsetpoint(NightSetPoint(), ts)
            @test isa(stp, Array{Float64})
            @test length(stp) == 96
        end
        @testset "Weather" begin
            for serie in [GTI, BHI, DHI, OutdoorTemperature]
                wt = loadweather(serie(), ts)
                @test isa(wt, Array{Float64})
                @test length(wt) == District.ntimesteps(ts)
            end
        end
    end


    @testset "Uncertainties" begin
        ts = TimeSpan(0, 1)

        @testset "Demands" begin
            demands = loadnoise(Demands(10, 1), ts)
            @test isa(demands, Array{Float64, 3})
            # test other timespan
            ts2 = TimeSpan(90, 6)
            demands = loadnoise(Demands(10, 1), ts2)
            @test isa(demands, Array{Float64, 3})
        end

        @testset "Fit & Forecast" begin
            for ξ in [Demands(10, 1), PVProduction(4, .15, 2, 0)]
                @test isa(District.elecload(ξ, 1), Union{Expr, Real})
                @test isa(District.nnoise(ξ), Int)

                # test fiting
                @test isa(District.fit(ξ, ts), WhiteNoise)

                # test generation of scenarios
                nscen = 10
                scen = District.genscenarios(ξ, ts, nscen)
                @test isa(scen, Array{Float64, 3})
                @test size(scen) == (District.ntimesteps(ts), nscen, District.nnoise(ξ))

                # test generation of forecast
                forecast = District.genforecast(ξ, ts)
                @test isa(forecast, Array{Float64, 2})
                @test size(forecast) == (District.ntimesteps(ts), District.nnoise(ξ))
            end
        end
    end


    devices = [Battery, ElecHotWaterTank, ThermalHotWaterTank, R6C2, R6C2, MicroCHP,
               ElecHeater, ThermalHeater]
    initid = ["bat0", "ehwt0", "twht0", "rt1988", "rt2012", "chp0", 6., 6.]
    @testset "Devices" begin
        for (Stock, ids) in zip(devices, initid)
            st = Stock(ids)
            @test isa(st, Stock)
            @test isa(District.elecload(st, 1), Union{Expr, Float64})
            @test isa(District.xbounds(st), Vector{Tuple{Float64, Float64}})
            @test isa(District.ubounds(st), Vector{Tuple{Float64, Float64}})
            @test length(District.xbounds(st)) == District.nstates(st)
            @test length(District.ubounds(st)) == District.ncontrols(st)
        end
    end


    @testset "Irradiation" begin
        env = R6C2("rt1988")
        ts = TimeSpan(0, 1)

        iwall, iwindow = get_irradiation(env, ts)

        for tab in [iwall, iwindow]
            @test isa(tab, Array{Float64})
            @test length(tab) == 96
        end
    end
end
