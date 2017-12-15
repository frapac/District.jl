
type Simulator
    ts::TimeSpan
    # SP Model
    model::StochDynamicProgramming.SPModel


    # configuration of assessment
    nassess::Int64

    generator::Symbol
end

function Simulator(model, idday, nassess, info;
                   generator=:real, structure=:pred, tol=.0, ndays=1)
    Simulator(Buildings.ID, model, idday, ndays, nassess, info,
              generator, structure, tol)
end


# TODO: import ndays in dataday
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
function simulate(model::StochDynamicProgramming.SPModel,
                  policy::Policy, x0::Vector{Float64},
                   scenario, forecast, real_dynamic, real_cost, predict::Function;
                   info=0, verbose=0)

    # Get number of timesteps:
    T = size(scenario, 1)
    nb_simulations = size(scenario, 2)

    stocks = zeros(T, nb_simulations, model.dimStates)
    controls = zeros(T, nb_simulations, model.dimControls)
    costs = zeros(nb_simulations)

    # Set first value of stocks equal to x0:
    for i in 1:nb_simulations
        stocks[1, i, :] = x0
    end

    p = Progress(T-1, 1)
    for t=1:T-1
        # If necessary, update forecast with last realizations:
        # update oracle:
        oracle = build_oracle(forecast)
        mpc_prob = mpcproblem(model, mpc, t, oracle, T)

        for k in 1:nb_simulations
            # get previous state:
            state_t = stocks[t, k, :]

            # find optimal control with MPC:
            if info == 0
                wt = collect(scenario[t, k, :])
                opt_control = mpccontrol(model, mpc_prob, state_t, wt)
            elseif info == 1
                wt = scenario[max(t-1, 1), k, :]
                opt_control = mpccontrol(model, mpc_prob, state_t, wt)
            elseif info == 3
                wt = predict(scenario, forecast, t, k)
                opt_control = mpccontrol(model, mpc_prob, state_t, wt)
            else
                opt_control = mpccontrol(model, mpc_prob, state_t)
            end

            alea_dh = vec(scenario[t, k, :])

            costs[k] += real_cost(t, state_t, opt_control, alea_dh)
            xf = real_dynamic(t, state_t, opt_control, alea_dh)

            stocks[t+1, k, :] = xf
            controls[t, k, :] = opt_control
        end
        next!(p)
    end

    for k = 1:nb_simulations
        costs[k] += realfinalcost(stocks[end, k, :])
    end
    return costs, stocks, controls, scenario
end
