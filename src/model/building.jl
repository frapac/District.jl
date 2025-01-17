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
const PENAL_TANK = 1.

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
    # links
    links::Vector{AbstractLink}
    # Inteface with graph
    conn::AbstractInterface
    # SP model
    model
end

House(ts::AbstractTimeSpan) = House(gensym(), ts, AbstractDevice[],
                                    AbstractUncertainty[], Billing(),
                                    (t, x, u, w) -> 0.,
                                    (t, x, u, w) -> 0.,
                                    AbstractLink[],
                                    NoneInterface(), nothing)
add!(h::House, dev::AbstractDevice) = push!(h.devices, dev)
add!(h::House, w::AbstractUncertainty) = push!(h.noises, w)
set!(h::House, p::AbstractPrice) = set!(h.billing, p)


nstocks(h::House) = sum(nstates.(h.devices))
ncontrols(h::House) = sum(ncontrols.(h.devices))
nnoises(h::House) = sum(nnoise.(h.noises))
ninjection(h::House) = 1
connectionsize(h::House) = ntimesteps(h.time) - 1

function build!(house::House, x0::Vector{Float64};
                info=:HD)
    buildload!(house)
    ntime = ntimesteps(house.time)
    xb = xbounds(house)
    ub = ubounds(house)

    reset!.(house.devices)
    buildlink!(house)

    # build probability laws
    laws = buildlaws(house)
    # build objective function
    costm = objective(house)
    # build final cost
    fcost = (info == :HD)? buildfcost(house) : final_cost_dh
    dynam = builddynamic(house)

    # build SP model
    spmodel = StochDynamicProgramming.LinearSPModel(ntime, ub,
                                                  x0, costm,
                                                  dynam,
                                                  tonoiselaws(laws),
                                                  info=info,
                                                  Vfinal=fcost)

    set_state_bounds(spmodel, xb)

    house.model = spmodel
end


# temp
function tonoiselaws(laws::WhiteNoise)
    NoiseLaw[NoiseLaw(laws[t].support', laws[t].probas) for t in 1:length(laws)]
end
function towhitenoise(laws::Vector{NoiseLaw})
    WhiteNoise([DiscreteLaw(w.support', w.proba) for w in laws])
end



################################################################################
# COST DEFINITION
################################################################################
function objective(house::House, xindex::Int=0, uindex::Int=0)
    bill = house.billing

    function costm(m, t, x, u, w)
        vals = JuMP.Variable[]
        coefs = Float64[]
        # add elec price
        if ~isa(bill.elec, NoneElecPrice)
            zel1 = @JuMP.variable(m, lowerbound=0)
            @constraint(m, zel1 >= bill.elec(t) * house.elecload(t, x, u, w))
            push!(vals, zel1)
            push!(coefs, 1.0)
            # add  injection price
            if ~isa(bill.injection, NoneElecPrice)
                @constraint(m, zel1 >= - bill.injection(t) * house.elecload(t, x, u, w))
            end
        end

        # add thermal comfort
        if ~isa(bill.comfort, NoneComfort)
            # get index of inner temperature
            itemp = getposition(house, R6C2) + xindex
            zth1 = @JuMP.variable(m, lowerbound=0)
            @constraint(m, zth1 >= -bill.comfort(t)*(x[itemp] - setpoint(bill.comfort, t) + 1))
            push!(vals, zth1)
            push!(coefs, 1.0)
        end

        # add gas price
        if ~isa(bill.gas, NoneGasPrice)
            zgas = @JuMP.variable(m, lowerbound=0)
            @constraint(m, zgas >= bill.gas(t)*house.gasload(t, x, u, w))
            push!(vals, zgas)
            push!(coefs, 1.0)
        end

        # TODO: move this part in dedicated function
        # add decomposition price
        if isa(house.conn, PriceInterface)
            # add < λ, F >
            u = m[:u]
            ifu = ncontrols(house) + uindex

            push!(vals, u[ifu])
            push!(coefs, house.conn.values[t])
            expr = JuMP.AffExpr(vals, coefs, 0.0)
            expr += JuMP.QuadExpr([u[ifu]], [u[ifu]], [1e-2], 0.)
        elseif isa(house.conn, FlowInterface)
            # add F = V^k
            u = m[:u]
            ifu = ncontrols(house) + uindex
            # TODO: do not consider these constraint in SDDP & assessment
            # add hard coupling constraint
            m.ext[:coupling] = @constraint(m, u[ifu] == house.conn.values[t])
            expr = JuMP.AffExpr(vals, coefs, 0.0)
            expr += JuMP.QuadExpr([u[ifu]], [u[ifu]], [1e-2], 0.)
        elseif isa(house.conn, QuadInterface)
            # add < λ, F > + norm(F - F^k)^2
            λk = house.conn.values[t]
            Fk = house.conn.penal[t]
            u = m[:u]
            ifu = ncontrols(house) + uindex
            push!(vals, u[ifu])
            push!(coefs, λk)
            expr = JuMP.AffExpr(vals, coefs, λk*Fk)
            expr += JuMP.QuadExpr([u[ifu]], [u[ifu]], [1e-2], 0.)
            # add quadratic penalty
            quadpenalty = u[ifu] + Fk
            expr += house.conn.τ/2. * quadpenalty^2
        else
            expr = JuMP.AffExpr(vals, coefs, 0.0)
        end

        return expr
    end
    return costm
end


################################################################################
# Final cost definition
################################################################################

# add a final penalty for considered devices
# WARNING: penalty should be given as an expression dependent of `x`.
# Ex:   :(2 * max(0, 2 - x))
function penalize!(house::House, dev::Type, p::Expr)
    pos = getposition(house, dev)
    p = MacroTools.postwalk(x -> x == :x ? Expr(:ref, :(xf), pos): x, p)
    add!(house.billing.finalcost, p)
end

# If final penalty is a max, we linearize it in optimization model
function addmax(model::JuMP.Model, z, ex::Expr)
    xf = model[:xf]
    for i in 2:endof(ex.args)
        f = eval(:(xf -> $(ex.args[i])))
        @constraint(model, z >= Base.invokelatest(f, xf))
    end
end

# WIP: definition of final cost for a house
# WARNING: TODO We should not use Base.invokelatest
function buildfcost(house::House)
    fcost = house.billing.finalcost
    function final_cost(model, m)
        alpha = m[:alpha]

        z = @JuMP.variable(m, [1:npenal(fcost)])
        # we add iteratively each max penalization as constraint
        for i in 1:npenal(fcost)
            addmax(m, z[i], fcost.ExprMax[i])
        end

        # add final cost as expression previously parsed,
        # with max linearized
        final = eval(:(z -> $(fcost.fcost_parsed)))
        @JuMP.constraint(m, alpha == Base.invokelatest(final, z))
    end
    return final_cost
end


# TODO: add DH final cost
function final_cost_dh(model, m)
    error("Unable to define final cost in DH: function currently broken")
    alpha = m[:alpha]
    # get number of random noises
    ns = model.noises[end-1].supportSize

    x = m[:x]
    u = m[:u]
    xf = m[:xf]
    @JuMP.variable(m, z1[1:ns] >= 0)
    @JuMP.constraint(m, z1[i=1:ns] .>= 2 - xf[2, i])
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
function parsebuilding(house::House, xindex::Int, uindex::Int, dt::Float64, params::Dict)
    exdyn = Expr[]
    for dev in house.devices
        dyn = parsedevice(dev, xindex, uindex, house.time.δt, params)
        xindex += nstates(dev)
        uindex += ncontrols(dev)

        push!(exdyn, dyn...)
    end
    return exdyn
end

function builddynamic(house::House)
    ntime = ntimesteps(house.time)
    params = Dict()
    if hasdevice(house, R6C2)
        pint, pext = get_irradiation(house)
        params["text"] = loadweather(OutdoorTemperature(), house.time)
        params["pint"] = pint
        params["pext"] = pext
    end

    xindex = 1
    uindex = 1
    exdyn = Expr(:vect)
    exdyn.args = parsebuilding(house, xindex, uindex, house.time.δt, params)

    return eval(:((t, x, u, w) -> $exdyn))
end



################################################################################
# LOAD DEFINITION
################################################################################
function buildload!(house::House, uindex::Int=1, windex::Int=1)
    ntime = ntimesteps(house.time)

    # build load corresponding to device
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
    for ξ in house.noises
        push!(excost.args, elecload(ξ, windex))
        windex += nnoise(ξ)
    end

    house.elecload = eval(:((t, x, u, w) -> $excost))
    house.gasload = eval(:((t, x, u, w) -> $exgas))
end


# get load of building only as Expression
function getload(house::House, uindex::Int=1, windex::Int=1)
    ntime = ntimesteps(house.time)

    # build load corresponding to device
    excost = Expr(:call, :+)
    for dev in house.devices
        load = elecload(dev, uindex)
        uindex += ncontrols(dev)
        push!(excost.args, load)
    end

    # build load corresponding to noise
    for ξ in house.noises
        push!(excost.args, elecload(ξ, windex))
        windex += nnoise(ξ)
    end
    return excost
end


################################################################################
# SIMULATION DEFINITION
################################################################################
"""Get real cost for simulation."""
function getrealcost(house::House)
    bill = house.billing

    # Note: we do not consider price decomposition cost in real cost
    # TODO: a lot of case are not considered here
    function real_cost(t, x, u, w)
        cost = 0.
        if ~isa(bill.elec, NoneElecPrice)
            flow  = house.elecload(t, x, u, w)
            cost += bill.elec(t)*max(0, flow)
        end
        if ~isa(bill.comfort, NoneComfort)
            itemp = getposition(house, R6C2)
            temp  = -x[itemp] + setpoint(bill.comfort, t) - 1
            cost += bill.comfort(t) * max(0, temp)
        end
        if ~isa(bill.gas, NoneGasPrice)
            cost += bill.gas(t)*house.gasload(t, x, u, w)
        end
        return cost
    end

    return real_cost
end
# TODO: final cost still hardcoded!!!
realfinalcost(xf) = PENAL_TANK*max(0., 5. - xf[2])

################################################################################
# LINKERS
################################################################################
join!(house::House, din::AbstractModel, dout::AbstractModel) = push!(house.links, Link(din, dout))

buildlink!(house::House, uindex::Int=0, windex::Int=0) = link!.(house, house.links, uindex, windex)

# TODO: use thermalload instead
# We can link EHWT only to DHW demands
function link!(house::House, hwt::ElecHotWaterTank, w::AbstractUncertainty, uind::Int, wind::Int)
    indw = windex(house, w) + wind
    push!(hwt.output.args, :(w[$indw]))
end

function link!(house::House, hwt::ThermalHotWaterTank, h::ThermalHeater, uind::Int, wind::Int)
    indu = uindex(house, h) + uind
    push!(hwt.output.args, :(u[$indu]))
end

function link!(house::House, hwt::ThermalHotWaterTank, chp::MicroCHP, uind::Int, wind::Int)
    indu = uindex(house, chp) + uind
    push!(hwt.input.args, :($(thermalload(chp, indu))))
end
function link!(house::House, hwt::ThermalHotWaterTank, burner::Burner, uind::Int, wind::Int)
    indu = uindex(house, burner) + uind
    push!(hwt.input.args, :($(thermalload(burner, indu))))
end

# add heater to R6C2
function link!(house::House, thm::R6C2, heat::AbstractHeater, uind::Int, wind::Int)
    indu = uindex(house, heat) + uind
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
getposition(house::House, dev::Type) = cumsum(nstates.(house.devices))[findfirst(isa.(house.devices, dev))]

"Speficy whether `house` has device with type `dev`."
hasdevice(house::House, dev::Type) = findfirst(isa.(house.devices, dev)) >= 1

"Speficy whether `house` has noise with type `dev`."
hasnoise(house::House, dev::Type) =  findfirst(isa.(house.noises, dev)) >= 1


################################################################################
# DECOMPOSITION
################################################################################
# TODO: here, setting Interface MUST happened just before building model
# TODO: raise error if devices added after Interface
function set!(house::House, conn::AbstractInterface)
    # if previous house connection was unset, add a linker to the building
    if isa(house.conn, NoneInterface)
        add!(house, conn.linker)
    else
        house.devices[end] = conn.linker
    end
    house.conn = conn
end

swap!(house::House, exch::Array{Float64}) = swap!(house.conn, exch)
flow!(house::House, flow::Array{Float64}) = flow!(house.conn, flow)


################################################################################
# PRINT
################################################################################
function show(io::IO, h::House)
    println("="^30)
    println("House $(h.name)")
    println("="^30)
    println("Available devices")
    println("-"^30)
    for d in h.devices
        println("- ", getname(d))
    end
    println("="^30)
end
