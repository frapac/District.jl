# Definition for random variables

export loadnoise, Demands

abstract type AbstractUncertainty end


function loadnoise end
function fit end


################################################################################
immutable Demands <: AbstractUncertainty end

function loadnoise(::Demands, ts::AbstractTimeSpan, idh=1)
    # select optimization scenarios subset
    nf, nt = 1, 1000
    ntime = ntimesteps(ts)
    ti, te = unravel(ts)

    db = getdb()

    elec = db["house$idh"]["De"]
    dhw  = db["house$idh"]["DHW"]

    noises = zeros(ntime, 1000, 2)

    # Select corresponding data :
    noises[:, :, 1] = elec[ti:te, nf:nt]
    noises[:, :, 2] = dhw[ti:te, nf:nt]
    return noises
end


################################################################################
immutable PVProduction <: AbstractUncertainty end

function loadnoise(::PVProduction, ts::AbstractTimeSpan) end
