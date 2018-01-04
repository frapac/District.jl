# Load deterministic data for District
# Deterministic data encompasses:
#   - prices (for electricity and comfort)
#   - setpoints (for internal temperatures)
#   - weather conditions (corresponding to year 2015)
#
# All timeseries are stored in a JLD database, stored in
#     data/weather.jld
# Other tariffs and septpoints are stored in
#     data/tariffs


export loadprice, EDFPrice, EPEXPrice
export loadsetpoint, NightSetPoint
export loadweather, OutdoorTemperature, GTI, DHI, BHI

abstract type AbstractData end

################################################################################
# Definition of prices
abstract type AbstractPrice <: AbstractData end
abstract type AbstractElecPrice <: AbstractPrice end

"""
    loadprice(p::AbstractPrice, ts::AbstractTimeSpan)

Load price values during date specified in TimeSpan `ts`. Return
price as `Vector{Float64}` with same size as `ts`.
"""
function loadprice end


# Off/On Peak tariffs
immutable EDFPrice <: AbstractElecPrice end
function loadprice(::EDFPrice, ts::AbstractTimeSpan)
    tariff = JSON.parsefile("data/tariffs/elec/edf.json")
    ntime = ntimesteps(ts)
    price = ones(ntime) * tariff["offpeak"]
    for nd=0:ts.ndays-1
        # TODO: ts is here hardcoded
        price[7*4+nd*96:23*4+nd*96] = tariff["onpeak"]
    end
    return price
end


# French EPEX price for day-ahead auction.
# Coming from:
# http://www.epexspot.com/en/
immutable EPEXPrice <: AbstractElecPrice end
loadprice(::EPEXPrice, ts::AbstractTimeSpan) = loaddata(ts, 16)


immutable Comfort <: AbstractPrice end
function loadprice(::Comfort, ts::AbstractTimeSpan)
    tariff = JSON.parsefile("data/tariffs/comfort/com0.json")
    return tariff["comfort"]
end


immutable EDFInjection <: AbstractElecPrice end
function loadprice(::EDFInjection, ts::AbstractTimeSpan)
    tariff = JSON.parsefile("data/tariffs/elec/edf.json")
    return tariff["inj"]
end

################################################################################
# Definition of temperature setpoints
abstract type AbstractSetPoint <: AbstractData end


"""
    loadsetpoint(p::AbstractSetPoint, ts::AbstractTimeSpan)

Load setpoint values during date specified in TimeSpan `ts`. Return
price as `Vector{Float64}` with same size as `ts`.
"""
function loadsetpoint end

# Night/Day setpoints
immutable NightSetPoint <: AbstractSetPoint end

function loadsetpoint(::NightSetPoint, ts::AbstractTimeSpan)
    setpoint = JSON.parsefile("data/tariffs/setpoints/night.json")
    ntime = ntimesteps(ts)
    price = ones(ntime) * setpoint["night"]
    for nd=0:ts.ndays-1
        # TODO: ts is here hardcoded
        price[6*4+nd*96:22*4+nd*96] = setpoint["day"]
    end
    return price
end



################################################################################
# Weather Data
abstract type AbstractWeatherData <: AbstractData end

"""
    loadweather(p::AbstractWeatherData, ts::AbstractTimeSpan)

Load weather values during date specified in TimeSpan `ts`. Return
price as `Vector{Float64}` with same size as `ts`.
"""
function loadweather end

immutable OutdoorTemperature <: AbstractWeatherData end
loadweather(::OutdoorTemperature, ts::AbstractTimeSpan) = loaddata(ts, 1) - 273.15

# Beam horizontal
immutable BHI <: AbstractWeatherData end
loadweather(::BHI, ts::AbstractTimeSpan) = loaddata(ts, 10)

# Diffuse horizontal
immutable DHI <: AbstractWeatherData end
loadweather(::DHI, ts::AbstractTimeSpan) = loaddata(ts, 11)

# Global tilted
immutable GTI <: AbstractWeatherData end
loadweather(::GTI, ts::AbstractTimeSpan) = loaddata(ts, 9)


################################################################################
# Generic
# load weather data from specified `column` in weather.jld
function loaddata(ts::AbstractTimeSpan, column::Int)
    weather = loadweather()
    ts, tf = unravel(ts)
    return weather[ts:tf, column]
end
