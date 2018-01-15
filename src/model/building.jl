################################################################################
# District.jl
################################################################################
# Generic models for  buildings.
# - Buildings currently implemented:
#    * House
# - Buildings are designed to lay on grid's nodes.
################################################################################
# TODO: add generic methods for AbstractBuilding

export House, add!, set!, link!, swap!
export nstocks, nnoises

# GENERIC TYPES
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
    billing::AbstractBilling
    # electricity load
    elecload::Function
    # gas load
    gasload::Function
    # Inteface with graph
    conn::AbstractInterface
    # SP model
    model
end

House(ts::AbstractTimeSpan) = House(gensym(), ts, AbstractDevice[],
                                    AbstractUncertainty[], Billing(),
                                    (t, x, u, w) -> 0.,
                                    (t, x, u, w) -> 0.,
                                    NoneInterface(), nothing)
add!(h::House, dev::AbstractDevice) = push!(h.devices, dev)
add!(h::House, w::AbstractUncertainty) = push!(h.noises, w)
set!(h::House, p::AbstractPrice) = set!(h.billing, p)


nstocks(h::House) = sum(nstates.(h.devices))
nnoises(h::House) = sum(nnoise.(h.noises))


function build!(house::House, x0::Vector{Float64})
    buildload!(house)
    ntime = ntimesteps(house.time)
    xb = xbounds(house)
    ub = ubounds(house)

    # build probability laws
    laws = buildlaws(house)
    # build objective function
    costm = objective(house)
    # build final cost
    fcost = final_cost
    dynam = builddynamic(house)

    # build SP model
    spmodel = StochDynamicProgramming.LinearSPModel(ntime, ub,
                                                  x0, costm,
                                                  dynam,
                                                  tonoiselaws(laws),
                                                  info=:HD,
                                                  Vfinal=fcost)

    set_state_bounds(spmodel, xb)
    # TODO: update model with connection

    house.model = spmodel
end


function updatemodel!(spmodel, conn::AbstractInterface)
    for t in spmodel.ntim
        updpb!(spmodel.solverInterface[t], conn, t, 1)
    end
end

# temp
function tonoiselaws(laws::WhiteNoise)
    NoiseLaw[NoiseLaw(laws[t].support', laws[t].probas) for t in 1:length(laws)]
end

################################################################################
# COST DEFINITION
################################################################################
# TODO: define more properly AffExpr
function objective(house::House)
    bill = house.billing

    function costm(m, t, x, u, w)
        vals = JuMP.Variable[]
        # add elec price
        if ~isa(bill.elec, NoneElecPrice)
            zel1 = @JuMP.variable(m, lowerbound=0)
            @constraint(m, zel1 >= bill.elec(t) * house.elecload(t, x, u, w))
            push!(vals, zel1)
            # add  injection price
            if ~isa(bill.injection, NoneElecPrice)
                @constraint(m, zel1 >= - bill.injection(t) * house.elecload(t, x, u, w))
            end
        end

        if ~isa(bill.comfort, NoneComfort)
            zth1 = @JuMP.variable(m, lowerbound=0)
            @constraint(m, zth1 >= -bill.comfort(t)*(x[4] - setpoint(bill.comfort, t) + 1))
            push!(vals, zth1)
        end

        if ~isa(bill.gas, NoneGasPrice)
            zgas = @JuMP.variable(m, lowerbound=0)
            @constraint(m, zgas >= bill.gas(t)*house.gasload(t, x, u, w))
            push!(vals, zgas)
        end

        coefs = ones(Float64, length(vals))
        return JuMP.AffExpr(vals, coefs, 0.0)
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

    return eval(:((t, x, u, w) -> $exdyn))
end



################################################################################
# LOAD DEFINITION
################################################################################
# TODO: build load on the fly
# TODO: clean hygiene of generated function
function buildload!(house::House)
    ntime = ntimesteps(house.time)

    # build load corresponding to device
    uindex = 1
    excost = Expr(:call, :+)
    exgas = Expr(:call, :+)
    for dev in house.devices
        load = elecload(dev, uindex)
        gas = gasload(dev, uindex)
        uindex += ncontrols(dev)
        push!(excost.args, load)
        push!(exgas.args, gas)
    end

    # build load corresponding to noise
    windex = 1
    for ξ in house.noises
        push!(excost.args, elecload(ξ, windex))
        windex += nnoise(ξ)
    end

    house.elecload = eval(:((t, x, u, w) -> $excost))
    house.gasload = eval(:((t, x, u, w) -> $exgas))
end



################################################################################
# SIMULATION DEFINITION
################################################################################
# TODO: clean definition of real cost
"""Get real cost for simulation."""
function getrealcost(house::House)
    bill = house.billing


    function real_cost(t, x, u, w)
        flow  = house.elecload(t, x, u, w)
        pelec = bill.elec(t)*max(0, flow)

        temp  = -x[4] + setpoint(bill.comfort, t) - 1
        pconfort = bill.comfort(t) * max(0, temp)

        return pelec + pconfort
    end

    return real_cost
end
realfinalcost(xf) = PENAL_TANK*max(0, 6 - xf[2])

################################################################################
# LINKERS
################################################################################
# TODO: use thermalload instead
# We can link EHWT only to DHW demands
function link!(house::House, hwt::ElecHotWaterTank, w::AbstractUncertainty)
    indw = windex(house, w)
    push!(hwt.output.args, :(w[$indw]))
end

function link!(house::House, hwt::ThermalHotWaterTank, h::ThermalHeater)
    indu = uindex(house, h)
    push!(hwt.output.args, :(u[$indu]))
end

function link!(house::House, hwt::ThermalHotWaterTank, chp::MicroCHP)
    indu = uindex(house, chp)
    push!(hwt.input.args, :($(thermalload(chp, indu))))
end

# add heater to R6C2
function link!(house::House, thm::R6C2, heat::AbstractHeater)
    indu = uindex(house, heat)
    push!(thm.input.args, :(u[$indu]))
end

################################################################################
# UTILS
################################################################################
get_irradiation(house::House) = get_irradiation(getdevice(house, R6C2), house.time)

# TODO: avoid side effect
"Get position of device in house.device."
getposition(house::House, d::AbstractDevice) = findfirst(house.devices, d)
getposition(house::House, w::AbstractUncertainty) = findfirst(house.noises, w)
xindex(house::House, d::AbstractDevice) = cumsum(nstates.(house.devices))[getposition(house, d)]
uindex(house::House, d::AbstractDevice) = cumsum(ncontrols.(house.devices))[getposition(house, d)]
windex(house::House, w::AbstractUncertainty) = cumsum(nnoise.(house.noises))[getposition(house, w)]

"Return first device with type `dev`."
getdevice(house::House, dev::Type) = house.devices[findfirst(isa.(house.devices, dev))]

"Speficy whether `house` has device with type `dev`."
hasdevice(house::House, dev::Type) = findfirst(isa.(house.devices, dev)) >= 1

"Speficy whether `house` has noise with type `dev`."
hasnoise(house::House, dev::Type) =  findfirst(isa.(house.noises, dev)) >= 1


################################################################################
# DECOMPOSITION
################################################################################
# TODO: here, setting Interface MUST happened just before building model
function set!(house::House, conn::PriceInterface)
    house.conn = conn
    add!(house, conn.linker)
    # TODO: update model
end

swap!(house::House, exch::Array{Float64}) = swap!(house.conn, exch)
