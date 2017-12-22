# Generate scenarios and forecasts for assessments

function genassessments(ts::AbstractTimeSpan, noises::Vector{AbstractUncertainty}, nscen::Int)
    # get noises dimensions
    nnoises = sum(nnoise.(noises))
    # allocate assessments
    assessment = zeros(Float64, ntimesteps(ts), nscen, nnoises)

    iw = 1
    for ξ in noises
        nw = nnoise(ξ) - 1
        assessment[:, :, iw:iw+nw] = genscenarios(ξ, ts, nscen)
        iw += nw + 1
    end
    return assessment
end



"""Simulate in-sample scenarios."""
function genscen(model::StochDynamicProgramming.SPModel, nscen::Int)
    StochDynamicProgramming.simulate_scenarios(model.noises, nscen)
end


################################################################################
# FORECASTING
################################################################################
function genforecast(ts::AbstractTimeSpan, noises::Vector{AbstractUncertainty})
    # get noises dimensions
    nnoises = sum(nnoise.(noises))
    # allocate assessments
    forecast = zeros(Float64, ntimesteps(ts), nnoises)

    iw = 1
    for ξ in noises
        nw = nnoise(ξ) - 1
        forecast[:, iw:iw+nw] = genforecast(ξ, ts)
        iw += nw + 1
    end
    return forecast
end
