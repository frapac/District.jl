################################################################################

function buildelechouse(ts::TimeSpan)
    house = House(ts)

    # Add devices
    devices = [Battery, ElecHotWaterTank, R6C2, ElecHeater]
    names = ["bat0", "ehwt0", "rt2012", 6.]
    dev = []
    for (Device, name) in zip(devices, names)
        d = Device(name)
        add!(house, d)
        push!(dev, d)
    end
    bat, hwt, thm, heat = dev
    # link heater to thermal envelope
    District.link!(house, thm, heat)

    # Add noises
    # import demands
    wdem = Demands(10, 1)
    add!(house, wdem)
    # import pv production
    wpv = PVProduction(1, .15, 20, 0)
    add!(house, wpv)

    # link hot water tank with hot water demand
    District.link!(house, hwt, wdem)

    # build objective: we penalize elec and thermal comfort.
    set!(house, EDFPrice(ts))
    set!(house, ComfortPrice(ts))

    # add initial pos


    return house
end

