
abstract type AbstractData end

################################################################################
# Definition of prices
abstract type AbstractPrice <: AbstractData end

function loadprice end

# Off/On Peak tariffs
immutable FixPrice <: AbstractPrice end


function loadprice(::FixPrice, iday, ndays)
    price = ones(ntime*ndays) * Params.P_ELEC
    for nd=0:ndays-1
        priceElec[7*4+nd*96:23*4+nd*96] = Params.P_ELEC_PEAK
    end
    return price
end


immutable EPEXPrice <: AbstractPrice end

function loadprice(::EPEXPrice, iday, ndays)
    #TODO: to implement
end




################################################################################
# Definition of temperature setpoints
abstract type AbstractSetPoint <: AbstractData end

immutable NightSetPoint <: AbstractSetPoint end

function loadsetpoint end


################################################################################
# Weather Data
abstract type AbstractWeatherData <: AbstractData end

immutable OutdoorTemperature <: AbstractWeatherData end
immutable Irradiation <: AbstractWeatherData end


