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

# TODO: add check consistency
export Battery, HotWaterTank, MicroCHP, R6C2

abstract type AbstractDevice end

"""
    elecload(d::AbstractDevice, uindex::Int)

Get impact of Device `d` on grid's load balance as an Expression.
`uindex` specified the starting index of devices control in grid's
dynamics.
"""
function elecload end

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

################################################################################
# Battery
struct Battery <: AbstractDevice
    name::String
    binf::Float64
    bup::Float64
    δb::Float64
    ρi::Float64
    ρe::Float64
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

elecload(bat::Battery, uindex) = :(u[$uindex] - u[$(uindex+1)])
nstates(bat::Battery) = 1
ncontrols(bat::Battery) = 2
xbounds(bat::Battery) = [(bat.binf, bat.bup)]
ubounds(bat::Battery) = [(0., bat.δb), (0., bat.δb)]


################################################################################
# EHWT
struct HotWaterTank <: AbstractDevice
    name::Symbol
    αt::Float64
    ηi::Float64
    ηe::Float64
    hmax::Float64
    power::Float64
end
function HotWaterTank(name::String)
    path = "$WD/data/devices/tank/$name.json"
    data = JSON.parsefile(path)

    HotWaterTank(name, data["ALPHA_H"], data["etain"], data["etaout"],
                 data["hmax"], data["power"])
end


# TODO: consistency with demands
function parsedevice(hwt::HotWaterTank, xindex::Int, uindex::Int, dt, p::Dict=Dict())
    dyn = [:($(hwt.αt)*x[$xindex] + $dt*($(hwt.ηi)*u[$uindex] - $(hwt.ηe)*w[2]))]
    return dyn
end

elecload(hwt::HotWaterTank, uindex::Int) = :(u[$uindex])

nstates(hwt::HotWaterTank) = 1
ncontrols(hwt::HotWaterTank) = 1
#TODO: dry bounds
xbounds(hwt::HotWaterTank) = [(0., hwt.hmax)]
ubounds(hwt::HotWaterTank) = [(0., hwt.power)]


################################################################################
# R6C2 model
struct R6C2 <: AbstractDevice
    name::Symbol
    altitude::Float64
    latitude::Float64

    Ssol::Float64
    Surf_window::Float64
    Surf_wall::Vector{Float64}
    H::Float64

    Ci::Float64
    Cw::Float64
    Giw::Float64
    Gie::Float64
    Gwe::Float64
    Ri::Float64
    Rs::Float64
    Rw::Float64
    Re::Float64

    Fv::Float64
    fframe::Float64
    albedo::Float64
    esp::Float64
    λe::Float64
    heater::Float64
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
    heater = params["heater"]
    R6C2(name, altitude, latitude,
         Ssol, Surf_window, Surf_wall, H,
         Ci, Cw, Giw, Gie, Gwe, Ri, Rs, Rw, Re,
         Fv, fframe, albedo, esp, le, heater)
end


function parsedevice(thm::R6C2, xindex::Int, uindex::Int, dt, p::Dict=Dict())
    dt2 = dt * 3600
    dyn = [:(
        x[$xindex]    + $dt2/$(thm.Cw)*($(thm.Giw)*(x[$(xindex+1)]-x[$xindex]) +
                                        $(thm.Gwe)*($(p["text"])[t]-x[$xindex]) +
                                        $(thm.Re)/($(thm.Re)+$(thm.Rw))*$(p["pext"])[t] +
                                        $(thm.Ri)/($(thm.Ri)+$(thm.Rs))*$(p["pint"])[t] +
                                        1000*$(thm.λe)*u[$xindex])), # wall's temperature
        :(x[$(xindex+1)] + $dt2/$(thm.Ci)*($(thm.Giw)*(x[$xindex]-x[$(xindex+1)]) +
                                         $(thm.Gie)*($(p["text"])[t]-x[$(xindex+1)]) +
                                         $(thm.Rs)/($(thm.Ri)+$(thm.Rs))*$(p["pint"])[t] +
                          1000*$(1-thm.λe)*u[$uindex]) # inner temperature
         )]
    return dyn
end

elecload(thm::R6C2, uindex::Int) = :(u[$uindex])
nstates(thm::R6C2) = 2
ncontrols(thm::R6C2) = 1
xbounds(thm::R6C2) = [(-50., 100.), (-50., 100.)]
ubounds(thm::R6C2) = [(0., thm.heater)]


################################################################################
# CHP model
# TODO: implement CHP
struct MicroCHP <: AbstractDevice
    name::String
    power::Float64
    yield::Float64
    eta_elec::Float64
    power_elec::Float64
    power_therm::Float64
    hwt::HotWaterTank
end
function MicroCHP(name)
    path = "$WD/data/devices/chp/$name.json"
    data = JSON.parsefile(path)
    power = data["CHP_POWER"]
    yield = data["CHP_YIELD"]
    eta = data["SHARE_ELEC"]
    hwt = HotWaterTank(data["watertank"])

    MicroCHP(name, power, yield, eta, power*yield*eta, power*yield*(1-eta), hwt)
end

function parsedevice(chp::MicroCHP, xindex::Int, uindex::Int, dt, p::Dict=Dict())
    dyn = [:($(chp.hwt.αt)*x[$xindex] + $dt*($(chp.hwt.ηi)*$(chp.power_therm)*u[$uindex] - $(chp.hwt.ηe)*(u[4]+w[2])))]
    return dyn
end

elecload(chp::MicroCHP, uindex::Int) = :(-u[$uindex]*$(chp.power_elec))
nstates(chp::MicroCHP) = 1
ncontrols(chp::MicroCHP) = 1
xbounds(chp::MicroCHP) = xbounds(chp.hwt)
ubounds(chp::MicroCHP) = [(0., 1.)]
