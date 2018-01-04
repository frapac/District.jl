# Definition for random variables

export loadnoise, Demands, PVProduction

abstract type AbstractUncertainty end


function loadnoise end
function nnoise end
function genscenarios end
function genforecast end
function fit end


################################################################################
immutable Demands <: AbstractUncertainty
    nbins::Int
    idhouse::Int
end

loadnoise(d::Demands, ts::AbstractTimeSpan) = _loadnoise(d, ts, 1, 1000)

function _loadnoise(d::Demands, ts::AbstractTimeSpan, nf, nt)
    nscen = nt - nf + 1
    idh = d.idhouse
    # select optimization scenarios subset
    ntime = ntimesteps(ts)
    wcycle = weekcycle(ts)

    db = getdb()

    elec = db["house$idh"]["De"]
    dhw  = db["house$idh"]["DHW"]

    noises = zeros(ntime, nscen, 2)

    # Select corresponding data :
    noises[:, :, 1] = elec[wcycle, nf:nt]
    noises[:, :, 2] = dhw[wcycle, nf:nt]
    return noises
end

elecload(d::Demands, windex) = :(w[$windex])
nnoise(d::Demands) = 2

function genscenarios(d::Demands, ts::AbstractTimeSpan, nscen::Int)
    @assert 1 <= nscen <= 1000
    return _loadnoise(d, ts, 1001, 1000 + nscen)
end

function genforecast(d::Demands, ts::AbstractTimeSpan)
    scen = genscenarios(d, ts, 1000)
    return mean(scen, 2)[:, 1, :]
end

fit(d::Demands, ts) = WhiteNoise(loadnoise(d, ts), d.nbins, KMeans())


################################################################################
immutable PVProduction <: AbstractUncertainty
    nbins::Int
    η::Float64
    area::Float64
    δ::Float64
end

#TODO: implement PV prod
#TODO: ensure tha PV prod is shared across all district

# quantize PV production
function loadnoise(::PVProduction, ts::AbstractTimeSpan) end

elecload(p::PVProduction, windex) = :(-w[$windex])
nnoise(p::PVProduction) = 1


function genscenarios(pv::PVProduction, ts::AbstractTimeSpan, nscen::Int)
    pvprod = production(pv, ts)

    scen = zeros(Float64, ntimesteps(ts), nscen, nnoise(pv))
    for n in 1:nscen
        scen[:, n, :] = genrandomproduction(pvprod, pv.δ)
    end
    return scen
end

"""Generate scenarios for PV with real laws."""
function genrandomproduction(val, δ)
    ntime = length(val)
    scen = zeros(ntime)

    for ii in 1:ntime
        ϵ = (ii-1)/(ntime-1)*δ*randn()
        scen[ii] = val[ii]*(1 + ϵ)
    end
    scen
end

genforecast(pv::PVProduction, ts::AbstractTimeSpan) = reshape(production(pv, ts), ntimesteps(ts), 1)
production(pv::PVProduction, ts::AbstractTimeSpan) = 0.001 * pv.η * pv.area * loadweather(GTI(), ts)


function fit(pv::PVProduction, ts::AbstractTimeSpan)
    mu = production(pv, ts)

    ntime = length(mu)
    laws = DiscreteLaw[]

    for t in 1:ntime
        # get variance
        ϵ = (t-1)/(ntime-1)*pv.δ
        # get optimal quantization of normal law
        quantization = normal_quantization(pv.nbins, 0, ϵ)
        # forecast is equal to μ (1 + N(0, ϵ))
        quantization[:, 2] = mu[t] * (1 + quantization[:, 2])
        push!(laws,  DiscreteLaw(quantization[:, 2], quantization[:, 1]))
    end

    return WhiteNoise(laws)
end


################################################################################
# UTILS
################################################################################
"""Get optimal quantization of normal law N(0, 1)."""
function normal_quantization(n)
    path = "data/quantization/"
    en = "_1_nopti"
    quantization = readdlm("$path$n$en")
    return quantization[1:end-1, 1:2]
end


"""Get optimal quantization of normal law N(μ, σ)"""
function normal_quantization(n, mu, sigma)
    quantization = normal_quantization(n)
    quantization[:, 2] = mu + sigma*quantization[:, 2]
    return quantization
end
