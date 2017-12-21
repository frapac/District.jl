# TODO: create SimulationResults
#

struct Simulator
    ts::TimeSpan
    # SP Model to simulate
    model::StochDynamicProgramming.SPModel
    # number of assessment scenarios
    scenarios
    realdynamic::Function
    realcost::Function
    realfinalcost::Function
end



function simulate(simulator::Simulator, policy::Policy)
    model = simulator.model

    # Get number of timesteps:
    ntime = size(scenario, 1)
    nsimu = size(scenario, 2)

    # Allocate
    stocks   = zeros(ntime, nsimu, model.dimStates)
    controls = zeros(ntime, nsimu, model.dimControls)
    costs    = zeros(nsimu)

    # Set first value of stocks equal to x0:
    for i in 1:nb_simulations
        stocks[1, i, :] = simulator.x0
    end

    @showprogress for t=1:ntime-1
        # update problem inside policy
        buildproblem!(policy, simulator.model, t)

        for k in 1:nsimu
            # get previous state:
            x = stocks[t, k, :]
            ξ = scenario[t, k, :]

            # compute decisions:
            u = policy(x, ξ)
            xf = simulator.realdynamic(t, x, u, ξ)

            costs[k] += simulator.realcost(t, x, u, ξ)
            stocks[t+1, k, :] = xf
            controls[t, k, :] = u
        end
    end
    for k = 1:nsimu
        costs[k] += simulator.realfinalcost(stocks[end, k, :])
    end
    return costs, stocks, controls
end
