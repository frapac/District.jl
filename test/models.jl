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

        @testset "Weather" begin
            for serie in [GTI, BHI, DHI, OutdoorTemperature]
                wt = loadweather(serie(), ts)
                @test isa(wt, Array{Float64})
                @test length(wt) == District.ntimesteps(ts)
            end
        end
    end

    @testset "Price" begin
        ts = TimeSpan(0, 1)
        @testset "Price" begin
            Prices = [EDFPrice, EPEXPrice, ComfortPrice, EngieGasPrice,
                      RecoursePrice, EDFInjection]
            for Tariff in Prices
                price = Tariff(ts)
                @test isa(price, District.AbstractPrice)
                @test isa(price(2), Float64)
            end
        end
        @testset "Setpoint" begin
            stp = NightSetPoint(ts)
            @test isa(stp, District.AbstractSetPoint)
        end
        @testset "Billing" begin
            bill = Billing()
            @test isa(bill, Billing)
            # test price loader
            set!(bill, EDFPrice(ts))
            set!(bill, EngieGasPrice(ts))
            set!(bill, EDFInjection(ts))
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

    @testset "Devices" begin
        # build params Dict to load devices
        params = Dict()
        params["text"] = zeros(Float64, 96)
        params["pint"] = zeros(Float64, 96)
        params["pext"] = zeros(Float64, 96)
        # add devices list
        devices = [Battery, ElecHotWaterTank, ThermalHotWaterTank, R6C2, R6C2, MicroCHP,
                ElecHeater, ThermalHeater, GraphConnection]
        initid = ["bat0", "ehwt0", "twht0", "rt1988", "rt2012", "chp0", 6., 6., 6.]
        # test devices
        for (Stock, ids) in zip(devices, initid)
            st = Stock(ids)
            @test isa(st, Stock)
            @test isa(District.elecload(st, 1), Union{Expr, Float64})
            @test isa(District.gasload(st, 1), Union{Expr, Float64})
            @test isa(District.thermalload(st, 1), Union{Expr, Float64})
            @test isa(District.xbounds(st), Vector{Tuple{Float64, Float64}})
            @test isa(District.ubounds(st), Vector{Tuple{Float64, Float64}})
            @test length(District.xbounds(st)) == District.nstates(st)
            @test length(District.ubounds(st)) == District.ncontrols(st)
            ex = District.parsedevice(st, 1, 1, .25, params)
            @test isa(ex, Vector{Expr})
        end
    end

    @testset "Graph Interface" begin
        conn = GraphConnection(6.)


        for Interface in [PriceInterface, FlowInterface]
            price = zeros(Float64, 96)
            p = Interface(price, conn)
            @test isa(p, District.AbstractInterface)
            price2 = ones(Float64, 96)
            District.swap!(p, price2)
            @test p.values == price2
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
