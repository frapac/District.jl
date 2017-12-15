# Implement deterministic data


export loadprice, EDFPrice, EPEXPrice
export loadsetpoint, NightSetPoint
export loadweather, OutdoorTemperature, GTI, DHI, BHI

abstract type AbstractData end

################################################################################
# Definition of prices
abstract type AbstractPrice <: AbstractData end

function loadprice end

# Off/On Peak tariffs
immutable EDFPrice <: AbstractPrice end
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


immutable EPEXPrice <: AbstractPrice end
loadprice(::EPEXPrice, ts::AbstractTimeSpan) = loaddata(ts, 16)




################################################################################
# Definition of temperature setpoints
abstract type AbstractSetPoint <: AbstractData end

function loadsetpoint end


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
function loadweather end

immutable OutdoorTemperature <: AbstractWeatherData end
loadweather(::OutdoorTemperature, ts::AbstractTimeSpan) = loaddata(ts, 1)

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
function loaddata(ts::AbstractTimeSpan, column::Int)
    weather = loadweather()
    ts, tf = unravel(ts)
    return weather[ts:tf, column]
end
