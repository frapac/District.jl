# Definition of generic devices

# TODO: load devices with JSON

abstract type AbstractDevice end


# Battery
struct Battery <: AbstractDevice
    name::Symbol
    binf::Float64
    bup::Float64
    δb::Float64
    ρi::Float64
    ρe::Float64
    αc::Float64
end


# EHWT
struct HotWaterTank <: AbstractDevice
    name::Symbol
    αt::Float64
    ηi::Float64
    ηe::Float64
end


# R6C2 model
struct R6C2 <: AbstractDevice
    Cw::Float64
    Giw::Float64
    Gwe::Float64
    Re::Float64
    Ri::Float64
    Rw::Float64
    Rs::Float64
    λe::Float64
end
