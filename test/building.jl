################################################################################
# District.jl
################################################################################
# Test District/models/buildings.jl
################################################################################

push!(LOAD_PATH, "..")

using Base.Test

using District


@testset "Building" begin
    @testset "Scratch" begin
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

        for d in [b, hwt, thm]
            idx = District.getposition(house, d)
            @test idx > 0
            @test isa(idx, Int)
            @test isa(District.xindex(house, d), Int)
            @test isa(District.uindex(house, d), Int)
        end
    end

    @testset "Sample buildings" begin
        ts = TimeSpan(0, 1)
        for HHouse in [ElecHouse, CHPHouse]
            house = load(ts, HHouse())
            @test isa(house, House)
        end
    end

    @testset "Graph interface" begin
        ts = TimeSpan(0, 1)
        # load house
        house = load(ts, ElecHouse())
        # define connection
        g = GraphConnection(6.)
        price = zeros(Float64, District.ntimesteps(ts))
        conn = PriceInterface(price, g)
        # add connection to house
        District.set!(house, conn)
        @test house.conn == conn
        @test District.hasdevice(house, GraphConnection)

        price2 = ones(Float64, District.ntimesteps(ts))
        District.swap!(house, price2)
        @test house.conn.price == price2
    end
end

