################################################################################
# District.jl
################################################################################
# Generic models for  buildings.
# - Buildings currently implemented:
#    * House
# - Buildings are designed to lay on grid's nodes.
################################################################################
# TODO: add generic methods for AbstractBuilding

export House, add!
export nstocks, nnoises

# GENERIC TYPES
# Definition of node in graph
abstract type AbstractNode end
abstract type AbstractBuilding <: AbstractNode end


################################################################################
# HOUSE
# TODO: penal tank hardcoded
PENAL_TANK = 1.

# TODO: add type for model. Shall we inherit from StochDynamicProgramming?
mutable struct House <: AbstractBuilding
    # name of House
    name::Symbol
    # time period considered
    time::AbstractTimeSpan
    # House's devices
    devices::Vector{AbstractDevice}
    # House's uncertainties
    noises::Vector{AbstractUncertainty}
    # House's prices
    prices::Vector{AbstractPrice}
    # SP model
    model
end

House(ts::AbstractTimeSpan) = House(gensym(), ts, AbstractDevice[],
                                    AbstractUncertainty[], AbstractPrice[], nothing)
add!(h::House, dev::AbstractDevice) = push!(h.devices, dev)
add!(h::House, w::AbstractUncertainty) = push!(h.noises, w)
add!(h::House, p::AbstractPrice) = push!(h.prices, p)


nstocks(h::House) = sum(nstates.(h.devices))
nnoises(h::House) = sum(nnoise.(h.noises))


function build!(house::House, x0::Vector{Float64})
    ntime = ntimesteps(house.time)
    xb = xbounds(house)
    ub = ubounds(house)


    laws = buildlaws(house)
    costm = objective(house)
    fcost = final_cost

    spmodel = StochDynamicProgramming.LinearSPModel(ntime, ub,
                                                  x0, costm,
                                                  dynam,
                                                  tonoiselaws(laws),
                                                  info=:HD,
                                                  Vfinal=fcost)

    set_state_bounds(spmodel, xb)
    house.model = spmodel
end


# temp
function tonoiselaws(laws::WhiteNoise)
    NoiseLaw[NoiseLaw(laws[t].support', laws[t].probas) for t in 1:length(laws)]
end

################################################################################
# COST DEFINITION
################################################################################
function objective(house::House)
    if hasdevice(house, R6C2)
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
    @JuMP.constraint(m, z1 >= 2. - xf[2])
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


################################################################################
# BOUNDS DEFINITION
################################################################################
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


################################################################################
# LAWS DEFINITION
################################################################################
function buildlaws(house::House)

    laws = WhiteNoise[]
    for ξ in house.noises
        push!(laws, fit(ξ, house.time))
    end

    # return product of laws
    return Scenarios.prodprocess(laws)
end


################################################################################
# DYNAMICS DEFINITION
################################################################################
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



################################################################################
# LOAD DEFINITION
################################################################################
function buildload(house::House)
    ntime = ntimesteps(house.time)

    # build load corresponding to device
    xindex = 1
    uindex = 1
    excost = Expr(:call, :+)
    for dev in house.devices
        load = elecload(dev, uindex)
        uindex += ncontrols(dev)
        push!(excost.args, load)
    end

    # build load corresponding to noise
    windex = 1
    for ξ in house.noises
        push!(excost.args, elecload(ξ, windex))
        windex += nnoise(ξ)
    end

    println(excost)
    @eval elecload(t, x, u, w) = $excost
    return elecload
end



################################################################################
# SIMULATION DEFINITION
################################################################################
# TODO: clean definition of real cost
"""Get real cost for simulation."""
function getrealcost(house::House)

    pel = loadprice(EDFPrice(), house.time)
    Tcons = loadsetpoint(NightSetPoint(), house.time)
    pin = loadprice(EDFInjection(), house.time)
    pth = loadprice(Comfort(), house.time)

    function real_cost(t, x, u, w)
        flow  = elecload(t, x, u, w)
        pelec = pel[t]*max(0, flow)

        temp  = -x[4] + Tcons[t] - 1
        pconfort = pth * max(0, temp)

        return pelec + pconfort
    end

    return real_cost
end
realfinalcost(xf) = PENAL_TANK*max(0, 6 - xf[2])


################################################################################
# UTILS
################################################################################
get_irradiation(house::House) = get_irradiation(getdevice(house, R6C2), house.time)

# TODO: avoid side effect
getdevice(house::House, dev::Type) = house.devices[findfirst(isa.(house.devices, dev))]
hasdevice(house::House, dev::Type) = findfirst(isa.(house.devices, dev)) >= 1
hasnoise(house::House, dev::Type) =  findfirst(isa.(house.noises, dev)) >= 1
