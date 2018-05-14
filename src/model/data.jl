################################################################################
# District.jl
################################################################################
# Load deterministic data for District.
# - Deterministic data encompasses:
#   * weather conditions (corresponding to year 2015)
# - All timeseries are stored in a JLD database, stored in
#     data/weather.jld
################################################################################


export loadweather, OutdoorTemperature, GTI, DHI, BHI

abstract type AbstractData end


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
# Convert Kelvin to Celsius Degree
loadweather(::OutdoorTemperature, ts::AbstractTimeSpan) = loaddata(ts, 1) - 273.15

immutable Humidity <: AbstractWeatherData end
loadweather(::Humidity, ts::AbstractTimeSpan) = loaddata(ts, 2)

### HORIZONTAL RADIATION ###
# Beam horizontal
immutable BHI <: AbstractWeatherData end
loadweather(::BHI, ts::AbstractTimeSpan) = loaddata(ts, 10)

# Diffuse horizontal
immutable DHI <: AbstractWeatherData end
loadweather(::DHI, ts::AbstractTimeSpan) = loaddata(ts, 11)

# Global horizontal
immutable GHI <: AbstractWeatherData end
loadweather(::GHI, ts::AbstractTimeSpan) = loaddata(ts, 12)

# Clear sky GHI
immutable ClearSkiGHI <: AbstractWeatherData end
loadweather(::ClearSkiGHI, ts::AbstractTimeSpan) = 4*loaddata(ts, 14)

### TILTED RADIATION ###
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
