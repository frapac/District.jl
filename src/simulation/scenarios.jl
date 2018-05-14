################################################################################
# District.jl
################################################################################
# Generate scenarios and forecasts for assessments.
################################################################################


export InSampleScenarios, OutSampleScenarios

################################################################################
# ASSESSMENT SCENARIOS
################################################################################
abstract type AbstractScenarioGenerator end

struct OutSampleScenarios <: AbstractScenarioGenerator
    shuffle::Bool
    seed::Int
end
OutSampleScenarios(;shuffle=true, seed=-1) = OutSampleScenarios(shuffle, seed)

struct InSampleScenarios <: AbstractScenarioGenerator
    seed::Int
end
InSampleScenarios(;seed=-1) = InSampleScenarios(seed)


"""
    genassessments(ts::AbstractTimeSpan, noises::Vector{AbstractUncertainty}, nscen::Int)

Generate `nscen` assessment scenarios for uncertainties `noises` over time period `ts`.

    genassessments(node::AbstractNode, nscen::Int)

Generate `nscen` assessment scenarios for uncertainties in Node `node`.

    genassessments(pb::Grid, nscen::Int)

Generate `nscen` assessment scenarios for uncertainties in Grid `pb`.
"""
function genassessments end

# TODO: add shuffling
function genassessments(ts::AbstractTimeSpan, noises::Vector{AbstractUncertainty}, nscen::Int, shuffle=false)
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
generate(gen::OutSampleScenarios, pb::AbstractGrid, nscen::Int) = genassessments(pb, nscen, shuff=gen.shuffle)
generate(gen::OutSampleScenarios, node::AbstractNode, nscen::Int) = genassessments(node, nscen)

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
genscen(n::AbstractNode, nscen) = genscen(n.model, nscen)

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

generate(gen::InSampleScenarios, pb::AbstractGrid, nscen::Int) = genscen(pb, nscen)
generate(gen::InSampleScenarios, node::AbstractNode, nscen::Int) = genscen(node, nscen)

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
