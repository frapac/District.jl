################################################################################
# District.jl
################################################################################
# Definition of generic devices.
# - Current implemented devices:
#       Battery, HotWaterTank, R6C2, MicroCHP
# - All Expression consider the notation of Optimal Control:
#    * `x` denotes the state
#    * `u` denotes the controls
#    * `w` denotes the uncertainties
################################################################################

export Battery, ElecHotWaterTank, MicroCHP, R6C2, ElecHeater, ThermalHotWaterTank,
       ThermalHeater
export GraphConnection

abstract type AbstractDevice <: AbstractModel end

"""
    elecload(d::AbstractDevice, uindex::Int)

Get impact of Device `d` on grid's load balance as an Expression.
`uindex` specified the starting index of devices control in grid's
dynamics.
"""
function elecload end

"""
    thermalload(d::AbstractDevice, uindex::Int)

Get thermal load outputed by device `d`.
`uindex` specified the starting index of devices control in grid's
dynamics.
"""
function thermalload end

"""
    gasload(d::AbstractDevice, uindex::Int)

Get gas load of device `d`.
"""
function gasload end

"""
    parsedevice(d::AbstractDevice, xindex::Int, uindex::Int, dt::Float64, p::Dict)

Get dynamics of Device `d` as an Expression.

# Arguments
* `d::AbstractDevice`: device to parse.
* `xindex::Int`: starting state index of device in grid.
* `uindex::Int`: starting control index of device in grid.
* `dt::Float64`: time delta of a timestep.
* `p::Dict`: optional parameters.
"""
function parsedevice end

"Get number of states in Device `d`."
function nstates end

"Get number of controls in Device `d`."
function ncontrols end

"Get states bounds for Device `d`."
function xbounds end

"Get controls bounds for Device `d`."
function ubounds end

"Specify whether device is stock"
function isstock end

"Get name of device."
function getname end

"""
    reset!(d::AbstractDevice)

Reset device."""
function reset! end

################################################################################
# Generic implementation
isstock(dev::AbstractDevice) = nstates(dev) > 0

# by default, elec, thermal and gas load are set to 0
elecload(d::AbstractDevice, uindex::Int) = :(0.)
thermalload(d::AbstractDevice, uindex::Int) = :(0.)
gasload(d::AbstractDevice, uindex::Int) = :(0.)
reset!(d::AbstractDevice) = nothing
getname(d::AbstractDevice) = "Undefined"

################################################################################
# Battery
struct Battery <: AbstractDevice
    name::String
    # battery lower bound
    binf::Float64
    # battery upper bound
    bup::Float64
    # max charge rate
    δb::Float64
    # charge yield
    ρi::Float64
    # discharge yield
    ρe::Float64
    # auto-discharge rate
    αc::Float64
end
function Battery(name::String)
    path = "$WD/data/devices/battery/$name.json"
    data = JSON.parsefile(path)

    Battery(name, data["BMIN"], data["BMAX"], data["DELTA_B_MAX"],
            data["rhoc"], data["rhod"], data["ALPHA_B"])
end


function parsedevice(bat::Battery, xindex::Int, uindex::Int, dt, p::Dict=Dict())
    dyn = [:($(bat.αc) * x[$xindex] + $dt*($(bat.ρi) * u[$uindex] - 1. / $(bat.ρe) * u[$(uindex+1)]))]
    return dyn
end

elecload(bat::Battery, uindex::Int) = :(u[$uindex] - u[$(uindex+1)])
nstates(bat::Battery) = 1
ncontrols(bat::Battery) = 2
xbounds(bat::Battery) = Tuple{Float64, Float64}[(bat.binf, bat.bup)]
ubounds(bat::Battery) = Tuple{Float64, Float64}[(0., bat.δb), (0., bat.δb)]
getname(bat::Battery) = "Battery"


################################################################################
# Electrical hot water tank
# EHWT has only a output and no input by construction.
struct ElecHotWaterTank <: AbstractDevice
    name::Symbol
    # auto-discharge rate
    αt::Float64
    # charge yield
    ηi::Float64
    # discharge yield
    ηe::Float64
    # tank max energy
    hmax::Float64
    # tank max power
    power::Float64
    # allowable temperature variation
    ΔT::Float64
    # output flow
    output::Expr
end
function ElecHotWaterTank(name::String)
    path = "$WD/data/devices/tank/$name.json"
    data = JSON.parsefile(path)

    dt = data["tempmax"] - data["tempmin"]
    hmax = 4.2*data["volume"]*dt / 3.6
    @assert hmax > 0
    # by default, output is null
    output = Expr(:call, +, 0)
    ElecHotWaterTank(name, data["ALPHA_H"], data["etain"], data["etaout"],
                 hmax, data["power"], dt, output)
end


function parsedevice(hwt::ElecHotWaterTank, xindex::Int, uindex::Int, dt, p::Dict=Dict())
    # convert l/h to W
    η = hwt.ηe * 4.2e-3 * 20
    dyn = [:($(hwt.αt)*x[$xindex] + $dt*($(hwt.ηi)*u[$uindex] - $(η)*$(hwt.output)))]
    return dyn
end

elecload(hwt::ElecHotWaterTank, uindex::Int) = :(u[$uindex])
thermalload(hwt::ElecHotWaterTank, uindex::Int) = :(0.)

reset!(hwt::ElecHotWaterTank)= (hwt.output.args = Any[:+])
nstates(hwt::ElecHotWaterTank) = 1
ncontrols(hwt::ElecHotWaterTank) = 1
xbounds(hwt::ElecHotWaterTank) = Tuple{Float64, Float64}[(0., hwt.hmax)]
ubounds(hwt::ElecHotWaterTank) = Tuple{Float64, Float64}[(0., hwt.power)]
getname(hwt::ElecHotWaterTank) = "EHWT"


################################################################################
# Thermal hot water tank
# THWT has only a output and no input by construction.
struct ThermalHotWaterTank <: AbstractDevice
    name::Symbol
    # auto-discharge rate
    αt::Float64
    # charge yield
    ηi::Float64
    # discharge yield
    ηe::Float64
    # tank max energy
    hmax::Float64
    # allowable temperature variation
    ΔT::Float64
    # input flow
    input::Expr
    # output flow
    output::Expr
end
function ThermalHotWaterTank(name::String)
    path = "$WD/data/devices/tank/$name.json"
    data = JSON.parsefile(path)

    dt = data["tempmax"] - data["tempmin"]
    hmax = 4.2*data["volume"]*dt / 3.6
    @assert hmax > 0
    # by default, output is null
    input = Expr(:call, +, 0)
    output = Expr(:call, +, 0)
    ThermalHotWaterTank(name, data["ALPHA_H"], data["etain"], data["etaout"],
                 hmax, dt, input, output)
end


function parsedevice(hwt::ThermalHotWaterTank, xindex::Int, uindex::Int, dt::Float64, p::Dict=Dict())
    dyn = [:($(hwt.αt)*x[$xindex] + $dt*($(hwt.ηi)*$(hwt.input) - $(hwt.ηe)*$(hwt.output)))]
    return dyn
end


function reset!(hwt::ThermalHotWaterTank)
    hwt.input.args = Any[:+]
    hwt.output.args = Any[:+]
end
nstates(hwt::ThermalHotWaterTank) = 1
ncontrols(hwt::ThermalHotWaterTank) = 0
xbounds(hwt::ThermalHotWaterTank) = Tuple{Float64, Float64}[(0., hwt.hmax)]
ubounds(hwt::ThermalHotWaterTank) = Tuple{Float64, Float64}[]
getname(hwt::ThermalHotWaterTank) = "THWT"

################################################################################
# R6C2 model
# by default, R6C2 has only input (corresponding to heating load)
struct R6C2 <: AbstractDevice
    name::Symbol
    altitude::Float64
    latitude::Float64

    # ground surface in m2
    Ssol::Float64
    # surface of window in m2
    Surf_window::Float64
    # surface of walls in m2
    Surf_wall::Vector{Float64}
    # house's height
    H::Float64

    # R6C2 parameters
    Ci::Float64
    Cw::Float64
    Giw::Float64
    Gie::Float64
    Gwe::Float64
    Ri::Float64
    Rs::Float64
    Rw::Float64
    Re::Float64

    # ventilation
    Fv::Float64
    fframe::Float64
    # wall's albedo
    albedo::Float64
    esp::Float64
    λe::Float64
    # input
    input::Expr
end
function R6C2(name)
    path = "$WD/data/devices/house/$name.json"
    params = JSON.parsefile(path)

    altitude = params["ALTITUDE"]
    latitude = params["LATITUDE"]

    # height
    H = params["H"]
    # Surface sol
    Sc = params["Sc"]
    N_floor = params["N_floor"]
    Ssol = Sc/N_floor
    # Surface windows
    Surf_window = sum(params["Surf_window"])
    # Surface walls
    Surf_wall = sqrt(Sc) * H * ones(5) - Surf_window
    Surf_wall[5] = Sc
    # isolation
    Utoit = params["Umur"]
    Uplancher = params["Umur"]
    Umur = params["Umur"]
    Ufen = params["Ufen"]

    ## R6C2 parameters:
    Ci = params["Ci"]
    Cw = params["Cw"]
    Gs = params["Gs"]
    Gi = 7.7*(2*Ssol + sqrt(Ssol)*H*4)
    Gw = Utoit*Ssol + Uplancher*Ssol + Umur*sqrt(Ssol)*H*4
    Ge = 25*(2*Ssol + sqrt(Ssol)*H*4)
    Gv = 150*1.2/3600*1000
    Rf = 1/(Ufen*Surf_window)
    Rinf = params["Rinf"]
    Gf = (Rinf + Rf)/(Rinf*Rf)
    Ri = 1/Gi
    Rs = 1/Gs
    Rw = 1/Gw
    Re = 1/Ge

    Giw = Gi*Gs/(Gi+Gs)
    Gie = Gv + Gf
    Gwe = Gw*Ge/(Gw+Ge)

    # lights
    Fv = params["Fv"]
    fframe = params["fframe"]
    albedo = params["albedo"]
    esp = params["esp"]
    le = params["lambdae"]
    # by default, input is null
    input = Expr(:call, +, 0)
    R6C2(name, altitude, latitude,
         Ssol, Surf_window, Surf_wall, H,
         Ci, Cw, Giw, Gie, Gwe, Ri, Rs, Rw, Re,
         Fv, fframe, albedo, esp, le, input)
end


function parsedevice(thm::R6C2, xindex::Int, uindex::Int, dt::Float64, p::Dict=Dict())
    dt2 = dt * 3600.
    # WARNING: here, the dynamic is express in Celsius degree instead of kWh
    dyn = [:(
        x[$xindex]    + $dt2/$(thm.Cw)*($(thm.Giw)*(x[$(xindex+1)]-x[$xindex]) +
                                        $(thm.Gwe)*($(p["text"])[t]-x[$xindex]) +
                                        $(thm.Re)/($(thm.Re)+$(thm.Rw))*$(p["pext"])[t] +
                                        $(thm.Ri)/($(thm.Ri)+$(thm.Rs))*$(p["pint"])[t] +
                                        1000*$(thm.λe)*$(thm.input))), # wall's temperature
        :(x[$(xindex+1)] + $dt2/$(thm.Ci)*($(thm.Giw)*(x[$xindex]-x[$(xindex+1)]) +
                                         $(thm.Gie)*($(p["text"])[t]-x[$(xindex+1)]) +
                                         $(thm.Rs)/($(thm.Ri)+$(thm.Rs))*$(p["pint"])[t] +
                                         1000*$(1-thm.λe)*$(thm.input)) # inner temperature
         )]
    return dyn
end

reset!(hwt::R6C2)= (hwt.input.args = Any[:+])
nstates(thm::R6C2) = 2
ncontrols(thm::R6C2) = 0
xbounds(thm::R6C2) = Tuple{Float64, Float64}[(-50., 100.), (-50., 100.)]
ubounds(thm::R6C2) = Tuple{Float64, Float64}[]
getname(thm::R6C2) = "Thermal Enveloppe"


################################################################################
# CHP and burners models
struct MicroCHP <: AbstractDevice
    name::String
    # max gas power
    power::Float64
    # thermal yield
    yield::Float64
    # proportion of elec
    eta_elec::Float64
    # max power elec
    power_elec::Float64
    # max thermal power
    power_therm::Float64
end
function MicroCHP(name)
    path = "$WD/data/devices/chp/$name.json"
    data = JSON.parsefile(path)
    power = data["CHP_POWER"]
    yield = data["CHP_YIELD"]
    eta = data["SHARE_ELEC"]
    MicroCHP(name, power, yield, eta, power*yield*eta, power*yield*(1-eta))
end

parsedevice(chp::MicroCHP, xindex::Int, uindex::Int, dt, p::Dict=Dict()) = Expr[]

elecload(chp::MicroCHP, uindex::Int) = :(-u[$uindex]*$(chp.power_elec))
thermalload(chp::MicroCHP, uindex::Int) = :(u[$uindex]*$(chp.power_therm))
gasload(chp::MicroCHP, uindex::Int) = :(u[$uindex] * $(chp.power))
nstates(chp::MicroCHP) = 0
ncontrols(chp::MicroCHP) = 1
xbounds(chp::MicroCHP) = Tuple{Float64, Float64}[]
ubounds(chp::MicroCHP) = Tuple{Float64, Float64}[(0., 1.)]
getname(chp::MicroCHP) = "μCHP"




################################################################################
# Heaters
abstract type AbstractHeater <: AbstractDevice end

struct ElecHeater <: AbstractHeater
    maxheating::Float64
end

parsedevice(h::ElecHeater, xindex::Int, uindex::Int, dt::Float64, p::Dict=Dict()) = Expr[]
elecload(h::ElecHeater, uindex::Int) = :(u[$uindex])
thermalload(h::ElecHeater, uindex::Int) = :(u[$uindex])
nstates(h::ElecHeater) = 0
ncontrols(h::ElecHeater) = 1
xbounds(h::ElecHeater) = Tuple{Float64, Float64}[]
ubounds(h::ElecHeater) = [(0., h.maxheating)]
getname(h::ElecHeater) = "Elec Heater"


# TODO: add dynamics of ThermalHeater
struct ThermalHeater <: AbstractHeater
    maxheating::Float64
end
parsedevice(h::ThermalHeater, xindex::Int, uindex::Int, dt::Float64, p::Dict=Dict()) = Expr[]
thermalload(h::ThermalHeater, uindex::Int) = :(u[$uindex])
nstates(h::ThermalHeater) = 0
ncontrols(h::ThermalHeater) = 1
xbounds(h::ThermalHeater) = Tuple{Float64, Float64}[]
ubounds(h::ThermalHeater) = [(0., h.maxheating)]
getname(h::ThermalHeater) = "Thermal Heater"


################################################################################
# Connection
abstract type AbstractConnection <: AbstractDevice end
struct GraphConnection <: AbstractConnection
    minkva::Float64
    maxkva::Float64
end

struct InOutGraphConnection <: GraphConnection

# By default, Node is prosumer
GraphConnection(kva::Float64) = GraphConnection(-kva, kva)

parsedevice(conn::GraphConnection, xindex::Int, uindex::Int, dt::Float64, p::Dict=Dict()) = Expr[]

elecload(conn::GraphConnection, uindex::Int) = :(-u[$uindex])
nstates(conn::GraphConnection) = 0
ncontrols(conn::GraphConnection) = 1
ncontrols(conn::InOutGraphConnection) = 2
xbounds(conn::GraphConnection) = Tuple{Float64, Float64}[]
ubounds(conn::GraphConnection) = Tuple{Float64, Float64}[(conn.minkva, conn.maxkva)]
getname(conn::GraphConnection) = "Import"


################################################################################
# Link between two devices
abstract type AbstractLink end
struct Link <: AbstractLink
    din::AbstractModel
    dout::AbstractModel
end

link!(n::AbstractNode, l::Link, uindex::Int, windex::Int) = link!(n, l.din, l.dout, uindex, windex)
