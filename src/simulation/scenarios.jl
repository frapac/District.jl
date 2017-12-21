

function genassessments(ts::AbstractTimeSpan, noises::Vector{AbstractUncertainty}, nscen::Int)
    # get noises dimensions
    nnoises = sum(nnoise.(noises))
    # allocate assessments
    assessment = zeros(Float64, ntimesteps(ts), nscen, nnoises)

    iw = 1
    for ξ in noises
        nw = nnoise(ξ)
        assessment[:, :, iw:iw+nw] = genscenarios(ξ, ts, nscen)
    end
    return assessment
end



"""Simulate in-sample scenarios."""
function genscen(model::StochDynamicProgramming.SPModel, nscen::Int)
    StochDynamicProgramming.simulate_scenarios(model.noises, nscen)
end


# TODO: clean
type DataDay
    scenario
    forecast
    dynamic
    function DataDay(simulator)
        idday = simulator.idday
        ndays = simulator.ndays
        nassess = simulator.nassess
        # Load demands and aleas:
        # Optim scenarios for forecast :
        demands_opt = get_aleas(96, scen="optim", idday=(idday-1)%7 + 1, ndays=ndays)
        # Assessment scenarios:
        demands = get_aleas(96, scen="assess", idday=(idday-1)%7 + 1, ndays=ndays)

        # Import weather condition:
        weather = import_data(1, ndays, start=(idday-1)*96+1)
        # Generate forecast with weather and optimization scenarios
        forecast = dayahead_forecast(simulator.model, demands_opt, weather, 1, 0,
                                     simulator.generator)
        # get dynamic corresponding to specified day
        realdynamic = get_dynamic(weather)

        scenario = genscen(simulator, weather, demands)

        return new(scenario, forecast, realdynamic)
    end
end
