# Definition of generic devices
using JSON

export Battery, HotWaterTank, CHP, R6C2
# TODO: load devices with JSON

abstract type AbstractDevice end

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

################################################################################
# EHWT
struct HotWaterTank <: AbstractDevice
    name::Symbol
    αt::Float64
    ηi::Float64
    ηe::Float64
end
function HotWaterTank(name::String)
    path = "$WD/data/devices/tank/$name.json"
    data = JSON.parsefile(path)

    HotWaterTank(name, data["ALPHA_H"], data["eta"], 1)
end

################################################################################
# R6C2 model
struct R6C2 <: AbstractDevice
    name::Symbol
    Ssol
    Surf_window
    Surf_wall
    H

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
end
function R6C2(name)
    path = "$WD/data/devices/house/$name.json"
    params = JSON.parsefile(path)

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
    R6C2(name, Ssol, Surf_window, Surf_wall, H,
         Ci, Cw, Giw, Gie, Gwe, Ri, Rs, Rw, Re,
         Fv, fframe, albedo, esp)
end



################################################################################
# CHP model
struct CHP <: AbstractDevice
    name
    power
    yield
    eta_elec
    power_elec
    power_therm
end
function CHP(name)
    path = "$WD/data/devices/chp/$name.json"
    data = JSON.parsefile(path)
    power = data["CHP_POWER"]
    yield = data["CHP_YIELD"]
    eta = data["SHARE_ELEC"]

    CHP(name, power, yield, eta, power*yield*eta, power*yield*(1-eta))
end
