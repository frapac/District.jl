################################################################################
# District.jl
################################################################################
# Load prices and billings for District.
# - Price data encompasses:
#   * prices (for electricity and comfort)
#   * setpoints (for internal temperatures)
# - Tariffs and septpoints are stored in
#     data/tariffs
################################################################################

################################################################################
export EDFPrice, EPEXPrice, ComfortPrice, EngieGasPrice
export NightSetPoint

# Definition of prices
abstract type AbstractPrice <: AbstractData end
abstract type AbstractSetPoint <: AbstractData end
abstract type AbstractBilling end
abstract type AbstractElecPrice <: AbstractPrice end
abstract type AbstractGasPrice <: AbstractPrice end
abstract type AbstractComfortPrice <: AbstractPrice end




################################################################################
# Elec Price
immutable NoneElecPrice <: AbstractElecPrice end

# Off/On Peak tariffs
struct EDFPrice <: AbstractElecPrice
    price::Vector{Float64}
end
function EDFPrice(ts::AbstractTimeSpan)
    tariff = JSON.parsefile("data/tariffs/elec/edf.json")
    ntime = ntimesteps(ts)
    price = ones(ntime) * tariff["offpeak"]
    for nd=0:ts.ndays-1
        # TODO: ts is here hardcoded
        price[7*4+nd*96:23*4+nd*96] = tariff["onpeak"]
    end
    return EDFPrice(price)
end
(p::EDFPrice)(t::Int) = p.price[t]


# French EPEX price for day-ahead auction.
# Coming from:
# http://www.epexspot.com/en/
struct EPEXPrice <: AbstractElecPrice
    price::Vector{Float64}
end
EPEXPrice(ts::AbstractTimeSpan) = EPEXPrice(loaddata(ts, 16))
(p::EPEXPrice)(t::Int) = p.price[t]


################################################################################
# Injection price
struct EDFInjection <: AbstractElecPrice
    price::Float64
end
function EDFInjection(ts::AbstractTimeSpan)
    tariff = JSON.parsefile("data/tariffs/elec/edf.json")
    return EDFInjection(tariff["inj"])
end
(p::EDFInjection)(t::Int) = p.price


################################################################################
# Gas price

immutable NoneGasPrice <: AbstractGasPrice end
struct EngieGasPrice <: AbstractGasPrice
    price::Float64
end
# TODO: dry
EngieGasPrice(ts::AbstractTimeSpan) = EngieGasPrice(0.06)
(p::EngieGasPrice)(t::Int) = p.price


################################################################################
# Setpoint

# Night/Day setpoints
struct NightSetPoint <: AbstractSetPoint
    setpoint::Vector{Float64}
end
function NightSetPoint(ts::AbstractTimeSpan)
    setpoint = JSON.parsefile("data/tariffs/setpoints/night.json")
    ntime = ntimesteps(ts)
    price = ones(ntime) * setpoint["night"]
    for nd=0:ts.ndays-1
        # TODO: ts is here hardcoded
        price[6*4+nd*96:22*4+nd*96] = setpoint["day"]
    end
    return NightSetPoint(price)
end


################################################################################
# Comfort Price
immutable NoneComfort <: AbstractComfortPrice end
struct ComfortPrice <: AbstractComfortPrice
    price::Float64
    setpoint::AbstractSetPoint
end
# TODO: setpoint is not necessarily NightSetPoint
function ComfortPrice(ts::AbstractTimeSpan)
    tariff = JSON.parsefile("data/tariffs/comfort/com0.json")
    ComfortPrice(tariff["comfort"], NightSetPoint(ts))
end
(p::ComfortPrice)(t::Int) = p.price
setpoint(c::ComfortPrice, t::Int) = c.setpoint.setpoint[t]


################################################################################
# Billing
mutable struct Billing <: AbstractBilling
    elec::AbstractElecPrice
    injection::AbstractElecPrice
    gas::AbstractGasPrice
    comfort::AbstractComfortPrice
end

Billing() = Billing(NoneElecPrice(), NoneElecPrice(), NoneGasPrice(), NoneComfort())
# type dispatch to set prices
set!(bill::Billing, p::AbstractElecPrice) = bill.elec = p
set!(bill::Billing, p::AbstractComfortPrice) = bill.comfort = p
set!(bill::Billing, p::AbstractGasPrice) = bill.gas = p
set!(bill::Billing, p::EDFInjection) = bill.injection = p
