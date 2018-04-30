################################################################################
# District.jl
################################################################################
# Implement classical housing profile.
# All profiles inherit from AbstractProfile
################################################################################


abstract type AbstractProfile end

export ElecHouse, CHPHouse, load
import JLD: load


"""
    load(ts::AbstractTimeSpan, prof::AbstractProfile)

Load house model corresponding to profile `prof`.
"""
function load end

function hasheater end


hasheater(p::AbstractProfile) = p.heater > 0


################################################################################
# Elec house
struct ElecHouse <: AbstractProfile
    # surface of PV [m2]
    pv::Float64
    # Heater [kW]
    heater::Float64
    # Battery profile
    bat::String
    # Thermal profile
    env::String
    # Demand profile
    idhouse::Int
    # Quantization
    nbinsdem::Int
end
ElecHouse(;pv=0, heat=3, bat="", env="rt2012", idhouse=1, nbins=10) = ElecHouse(pv, heat, bat, env, idhouse, nbins)


function load(ts::TimeSpan, prof::ElecHouse)
    house = House(ts)

    # Add battery
    if prof.bat != ""
        d = Battery(prof.bat)
        add!(house, d)
    end

    hwt = ElecHotWaterTank("ehwt0")
    add!(house, hwt)
    penalize!(house, ElecHotWaterTank, :(1. * max(0, 2 -x)))

    if hasheater(prof)
        thm = R6C2(prof.env)
        add!(house, thm)

        heat = ElecHeater(prof.heater)
        add!(house, heat)

        # link heater to thermal envelope
        join!(house, thm, heat)
    end

    # import demands
    wdem = Demands(prof.nbinsdem, prof.idhouse)
    add!(house, wdem)

    # import pv production
    if prof.pv > 0.
        wpv = PVProduction(1, .15, prof.pv, 0)
        add!(house, wpv)
    end

    # link hot water tank with hot water demand
    join!(house, hwt, wdem)

    # build objective: we penalize elec and thermal comfort.
    set!(house, EDFPrice(ts))

    if hasheater(prof)
        set!(house, ComfortPrice(ts))
    end

    return house
end

################################################################################
# CHP house
struct CHPHouse <: AbstractProfile
    chp::String
    # Heater [kW]
    heater::Float64
    # Battery profile
    bat::String
    # Thermal profile
    env::String
    # Demand profile
    idhouse::Int
    # Quantization
    nbinsdem::Int
end
CHPHouse(;chp="chp0", heat=6., bat="", env="rt1988", idhouse=1, nbins=10) = CHPHouse(chp, heat, bat, env, idhouse, nbins)


function load(ts::TimeSpan, prof::CHPHouse)
    house = House(ts)

    chp = MicroCHP(prof.chp)
    add!(house, chp)
    # Add battery
    if prof.bat != ""
        d = Battery(prof.bat)
        add!(house, d)
    end

    hwt = ThermalHotWaterTank("twht0")
    add!(house, hwt)
    penalize!(house, ThermalHotWaterTank, :(1. * max(0, 5 -x)))

    thm = R6C2(prof.env)
    add!(house, thm)

    heat = ThermalHeater(prof.heater)
    add!(house, heat)

    boiler = Burner(4., yield=.9)
    add!(house, boiler)

    join!(house, thm, heat)
    join!(house, hwt, heat)
    join!(house, hwt, chp)
    join!(house, hwt, boiler)

    # import demands
    wdem = Demands(prof.nbinsdem, prof.idhouse)
    add!(house, wdem)

    set!(house, EDFPrice(ts))
    set!(house, ComfortPrice(ts))
    set!(house, District.EngieGasPrice(ts))

    return house
end
