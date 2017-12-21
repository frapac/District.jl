# Definition for random variables

# TODO: add quantization bins in Uncertainty?

export loadnoise, Demands

abstract type AbstractUncertainty end


function loadnoise end


################################################################################
immutable Demands <: AbstractUncertainty end

loadnoise(d::Demands, ts::AbstractTimeSpan, idh=1) = _loadnoise(d, ts, idh, 1, 1000)

function _loadnoise(::Demands, ts::AbstractTimeSpan, idh, nf, nt)
    # select optimization scenarios subset
    ntime = ntimesteps(ts)
    wcycle = weekcycle(ts)

    db = getdb()

    elec = db["house$idh"]["De"]
    dhw  = db["house$idh"]["DHW"]

    noises = zeros(ntime, 1000, 2)

    # Select corresponding data :
    noises[:, :, 1] = elec[wcycle, nf:nt]
    noises[:, :, 2] = dhw[wcycle, nf:nt]
    return noises
end

elecload(d::Demands, windex) = :(w[$windex])
nnoise(d::Demands) = 2

function genscenarios(d::Demands, ts::AbstractTimeSpan, nscen::Int, idh=1)
    @assert nscen <= 1000
    return _loadnoise(d, ts, idh, 1001, 1000 + nscen)
end


################################################################################
immutable PVProduction <: AbstractUncertainty end

#TODO: implement PV prod

function loadnoise(::PVProduction, ts::AbstractTimeSpan) end
elecload(p::PVProduction, windex) = :(-w[$windex])
nnoise(p::PVProduction) = 1

#TODO: adapt
"""Generate scenarios for PV with real laws."""
function _genperturbations(val, δ)
    error("Deprecated")
    ntime = length(val)
    scen = zeros(ntime)

    for ii in 1:ntime
        ϵ = (ii-1)/(ntime-1)*δ*randn()
        scen[ii] = val[ii]*(1 + ϵ)
    end
    scen
end
