
# TODO: dry P_TH
# TODO: define p_inj
export House, add!
export nstocks

abstract type AbstractBuilding end


struct House <: AbstractBuilding
    name::Symbol
    time::AbstractTimeSpan
    devices::Vector{AbstractDevice}
    noises::Vector{AbstractUncertainty}
    prices::Vector{AbstractPrice}
    model
end

House(ts::AbstractTimeSpan) = House(gensym(), ts, AbstractDevice[],
                                    AbstractUncertainty[], AbstractPrice[], nothing)
add!(h::House, dev::AbstractDevice) = push!(h.devices, dev)
add!(h::House, w::AbstractUncertainty) = push!(h.noises, w)
add!(h::House, p::AbstractPrice) = push!(h.prices, p)


nstocks(h::House) = sum(nstates.(h.devices))


function objective(house::House)
    if has(house, R6C2)
        return _objectivethermal(house)
    else
        return _objectiveelec(house)
    end
end

# TODO: clean objective dispatch
# FIXME: must be run after buildload
function _objectiveprice(house::House)
    p_elec = loadprice(EDFPrice(), ts)
    p_inj = loadprice(EDFInjection(), ts)


    function costm(m, t, x, u, w)
        zel1 = @JuMP.variable(m, lowerbound=0)
        @constraint(m, zel1 >= p_elec[t] * elecload(t, x, u, w))
        @constraint(m, zel1 >= - p_inj * elecload(t, x, u, w))

        return JuMP.AffExpr(zel1)
    end
    return costm
end

function _objectivethermal(house::House)
    p_elec = loadprice(EDFPrice(), house.time)
    p_therm = loadsetpoint(NightSetPoint(), house.time)
    p_inj = loadprice(EDFInjection(), house.time)
    p_th = loadprice(Comfort(), house.time)


    function costm(m, t, x, u, w)
        zel1 = @JuMP.variable(m, lowerbound=0)
        @constraint(m, zel1 >= p_elec[t] * elecload(t, x, u, w))
        @constraint(m, zel1 >= - p_inj * elecload(t, x, u, w))

        # TODO: handle HVAC
        zth1 = @JuMP.variable(m, lowerbound=0)
        @constraint(m, zth1 >= -p_th*(x[4] - p_therm[t] + 1))

        return JuMP.AffExpr(zel1 + zth1)
    end
    return costm
end


function final_cost(model, m)
    alpha = m[:alpha]
    #= w = JuMP.getvariable(m, :w) =#
    x = m[:x]
    u = m[:u]
    xf = m[:xf]
    @variable(m, cost)
    z1 = @JuMP.variable(m, lowerbound=0)
    @JuMP.constraint(m, z1 >= 6 - xf[2])
    @JuMP.constraint(m, alpha == PENAL_TANK*z1)
end


function final_cost_dh(model, m)
    alpha = m[:alpha]
    # get number of random noises
    ns = model.noises[end-1].supportSize

    x = m[:x]
    u = m[:u]
    xf = m[:xf]
    @JuMP.variable(m, z1[1:ns] >= 0)
    @JuMP.constraint(m, z1[i=1:ns] .>= 6 - xf[2, i])
    for i in 1:ns
        @JuMP.constraint(m, alpha[i] == PENAL_TANK*z1[i])
    end
end

#TODO
function xbounds(house::House)
    xb= Tuple{Float64, Float64}[]
    for dev in house.devices
        xb = vcat(xb, xbounds(dev))
    end
    return xb
end

function ubounds(house::House)
    ub= Tuple{Float64, Float64}[]
    for dev in house.devices
        ub = vcat(ub, ubounds(dev))
    end
    return ub
end


function buildlaw(house::House)
    demands = loadnoise(Demands(), house.time)
    laws = WhiteNoise(demands, nbins, KMeans())
    # TODO
    #= lawpv = #TODO =#
    #= laws = prodlaw(lawdemands, lawpv) =#
end


function builddynamic(house::House)
    x0=[.6, 6, 16, 16]
    ntime = ntimesteps(house.time)
    pint, pext = get_irradiation(house)
    params = Dict()
    params["text"] = loadweather(OutdoorTemperature(), house.time)
    params["pint"] = pint
    params["pext"] = pext

    xindex = 1
    uindex = 1
    exdyn = Expr(:vect)
    for dev in house.devices
        dyn = parsedevice(dev, xindex, uindex, house.time.δt, params)
        xindex += nstates(dev)
        uindex += ncontrols(dev)

        for d in dyn
            push!(exdyn.args, d)
        end
    end

    @eval dynam(t, x, u, w) = $exdyn
    return dynam
end

function buildload(house::House)
    ntime = ntimesteps(house.time)

    xindex = 1
    uindex = 1
    excost = Expr(:call, :+)
    for dev in house.devices
        load = elecload(dev, uindex)
        uindex += ncontrols(dev)
        push!(excost.args, load)
    end

    @eval elecload(t, x, u, w) = $excost
    return elecload
end


function dynamic(t, x, u, w)
    return [Params.ALPHA_B * x[1] + Params.DT*(Params.ρc*u[1] - 1/Params.ρd*u[2]), # Battery
            x[2] + Params.DT*(u[4]*Params.ETA_R - w[2]), # Hot water tank
            x[3] + Params.DT*3600/Params.Cw*(Params.Giw*(x[4]-x[3]) +
                                        Params.Gwe*(w[3]-x[3]) +
                                        Params.Re/(Params.Re+Params.Rw)*w[4] +
                                        Params.Ri/(Params.Ri+Params.Rs)*α*w[4] +
                                        1000*.35*u[3]), # wall's temperature
            x[4] + Params.DT*3600/Params.Ci*(Params.Giw*(x[3]-x[4]) +
                                        Params.Gie*(w[3]-x[4]) +
                                        Params.Rs/(Params.Ri+Params.Rs)*α*w[4] +
                                        1000*.65*u[3])] # inner temperature
end




"""Get real cost for assessment."""
function get_real_cost(day)
    priceElec, Tcons = get_reference()

    function real_cost(t, x, u, w)
        flow  = (w[1] + u[4]- w[5] - u[2] + u[1] + u[3])
        pelec = priceElec[t]*max(0, flow)

        temp  = -x[4] + Tcons[t] - 1
        pconfort = Params.P_TH * max(0, temp)

        return pelec + pconfort
    end

    return real_cost
end
realfinalcost(xf) = PENAL_TANK*max(0, 6 - xf[2])



# UTILS
get_irradiation(house::House) = get_irradiation(getdevice(house, R6C2), house.time)

# TODO: avoid side effect
getdevice(house::House, dev::Type) = house.devices[findfirst(isa.(house.devices, dev))]
has(house::House, dev::Type) = findfirst(isa.(house.devices, dev)) >= 1
