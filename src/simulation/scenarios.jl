
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

function genscen(simulator, weather, demands)
    if simulator.generator == :real
        return _genscen(weather, demands, simulator.nassess,
                        simulator.model.dimNoises, ndays=simulator.ndays)
    elseif simulator.generator == :random
        return _genscen(simulator.model, simulator.nassess)
    elseif simulator.generator == :pv
        return _genscen(simulator, weather, simulator.nassess)
    end
end


"""Simulate demands scenarios."""
function _genscen(weather::AleasBuilding,
                       demands::Array{Float64, 3},
                       nscen::Int, naleas::Int; ndays=1)
    step_per_day = 96*ndays
    scenarios = zeros(step_per_day, nscen, naleas)
    scenarios[:, :, 1:naleas-3] = demands[:, 1:nscen, :]
    for i in 1:nscen
        scenarios[:, i, naleas-2] = weather.outdoor_temperature[1:96*ndays]
        scenarios[:, i, naleas-1] = weather.radiative_external_gain[1:96*ndays]
        scenarios[:, i, naleas] = weather.prod_pv[1:96*ndays]
    end
    return scenarios
end


"""Simulate in-sample scenarios."""
function _genscen(model::StochDynamicProgramming.SPModel, nscen::Int)
    StochDynamicProgramming.simulate_scenarios(model.noises, nscen)
end


"""Simulate PV scenarios."""
function _genscen(simulator, weather, nscen)
    scenarios = _genscen(simulator.model, nscen)
    # FIXME: we assume that PV values are at the end of the array
    for ii in 1:nscen
        scenarios[:, ii, end] = _genperturbations(weather.prod_pv, simulator.tol_pv)
    end
    return scenarios
end

"""Generate scenarios for PV with real laws."""
function _genperturbations(val, δ)
    ntime = length(val)
    scen = zeros(ntime)

    for ii in 1:ntime
        ϵ = (ii-1)/(ntime-1)*δ*randn()
        scen[ii] = val[ii]*(1 + ϵ)
    end
    scen
end
