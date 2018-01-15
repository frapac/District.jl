################################################################################
# District.jl
################################################################################
# Implement classical housing profile.
# All profiles inherit from AbstractProfile
################################################################################


abstract type AbstractProfile end

import JLD: load


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
end
ElecHouse(;pv=0, heat=3, bat="", env="rt2012", idhouse=1) = ElecHouse(pv, heat, bat, env, idhouse)


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
    District.link!(house, thm, heat)

    # import demands
    wdem = Demands(10, prof.idhouse)
    add!(house, wdem)

    # import pv production
    if prof.pv > 0.
        wpv = PVProduction(1, .15, prof.pv, 0)
        add!(house, wpv)
    end

    # link hot water tank with hot water demand
    District.link!(house, hwt, wdem)

    # build objective: we penalize elec and thermal comfort.
    set!(house, EDFPrice(ts))
    set!(house, ComfortPrice(ts))

    return house
end

