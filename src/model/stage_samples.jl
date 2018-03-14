################################################################################
# District.jl
################################################################################
# Implement classical housing profile.
# All profiles inherit from AbstractProfile
################################################################################


abstract type AbstractProfile end

export ElecHouse, CHPHouse, load, loadMultipleHouse
import JLD: load

srand(2713)

"""
    load(ts::AbstractTimeSpan, prof::AbstractProfile)

Load house model corresponding to profile `prof`.
"""
function load end


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

function loadMultipleHouse(ts::TimeSpan, nnodes::Int64)
    houseArray = Array{AbstractNode,1}(nnodes)
    for i in 1:nnodes
        solarpv = bitrand(1)
        if solarpv[1]
            houseArray[i] = load(ts, ElecHouse(pv=rand(1:10), heat=rand(2:6), bat="bat0"^rand(0:1), nbins=1))
        else 
            houseArray[i] = load(ts, ElecHouse(pv=0, heat=rand(2:6), bat="", nbins=1))
        end
    end

    return houseArray
end


function load(ts::TimeSpan, prof::ElecHouse)
    house = House(ts)

    # Add battery
    if prof.bat != ""
        d = Battery(prof.bat)
        add!(house, d)
    end

    hwt = ElecHotWaterTank("ehwt0")
    add!(house, hwt)

    thm = R6C2(prof.env)
    add!(house, thm)

    heat = ElecHeater(prof.heater)
    add!(house, heat)

    # link heater to thermal envelope
    join!(house, thm, heat)

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
    set!(house, ComfortPrice(ts))

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
end
CHPHouse(;chp="chp0", heat=6., bat="", env="rt1988", idhouse=1) = CHPHouse(chp, heat, bat, env, idhouse)


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

    thm = R6C2(prof.env)
    add!(house, thm)

    heat = ThermalHeater(prof.heater)
    add!(house, heat)

    join!(house, thm, heat)
    join!(house, hwt, heat)
    join!(house, hwt, chp)

    # import demands
    wdem = Demands(10, prof.idhouse)
    add!(house, wdem)

    set!(house, EDFPrice(ts))
    set!(house, ComfortPrice(ts))
    set!(house, District.EngieGasPrice(ts))

    return house
end
