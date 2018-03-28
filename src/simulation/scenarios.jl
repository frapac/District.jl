################################################################################
# District.jl
################################################################################
# Generate scenarios and forecasts for assessments.
################################################################################


################################################################################
# ASSESSMENT SCENARIOS
################################################################################

"""
    genassessments(ts::AbstractTimeSpan, noises::Vector{AbstractUncertainty}, nscen::Int)

Generate `nscen` assessment scenarios for uncertainties `noises` over time period `ts`.

    genassessments(node::AbstractNode, nscen::Int)

Generate `nscen` assessment scenarios for uncertainties in Node `node`.

    genassessments(pb::Grid, nscen::Int)

Generate `nscen` assessment scenarios for uncertainties in Grid `pb`.
"""
function genassessments end

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

# overload to generate assessment scenarios directly with node
genassessments(node::AbstractNode, nscen::Int) = genassessments(node.time, node.noises, nscen)

function genassessments(pb::Grid, nscen::Int; shuff=true)
    # get total number of uncertainties
    nw = sum(nnoises.(pb.nodes))

    scenarios = zeros(Float64, ntimes(pb), nscen, nw)

    iw = 1

    for node in pb.nodes
        nwnode = nnoises(node)
        nodescen = genassessments(node, nscen)
        if shuff
            nodescen[:] = nodescen[:, randperm(nscen), :]
        end
        scenarios[:, :, iw:iw+nwnode-1] = nodescen
        iw += nwnode
    end
    return scenarios
end


"""Simulate in-sample scenarios."""
function genscen(model::StochDynamicProgramming.SPModel, nscen::Int)
    StochDynamicProgramming.simulate_scenarios(model.noises, nscen)
end
function genscen(pb::Grid, nscen::Int)
    # get total number of uncertainties
    nw = sum(nnoises.(pb.nodes))
    scenarios = zeros(Float64, ntimes(pb), nscen, nw)

    iw = 1
    for node in pb.nodes
        nwnode = nnoises(node)
        scenarios[:, :, iw:iw+nwnode-1] = genscen(node.model, nscen)
        iw += nwnode
    end
    return scenarios

end


################################################################################
# FORECASTING
################################################################################
"""
    genforecast(ts::AbstractTimeSpan, noises::Vector{AbstractUncertainty})

Generate a forecast for uncertainties `noises` over time period `ts`.
"""
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
