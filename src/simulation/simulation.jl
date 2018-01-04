# Run simulation on given model and scenarios

export Simulator, SimulationResult
import Base: show


################################################################################
# Simulation results
struct SimulationResult
    costs::Vector{Float64}
    stocks::Array{Float64, 3}
    controls::Array{Float64, 3}
end

# TODO: .25 hardcoded
function show(io::IO, res::SimulationResult)
    @printf("Costs: %.4f ± %.4f €\n", mean(.25res.costs), std(.25res.costs))
end

################################################################################
# Simulator
struct Simulator
    ts::TimeSpan
    # SP Model to simulate
    model::StochDynamicProgramming.SPModel
    # number of assessment scenarios
    scenarios::Array{Float64, 3}
    realdynamic::Function
    realcost::Function
    realfinalcost::Function
end
function Simulator(n::AbstractNode, nassess::Int)
    ts = n.time
    scen = genassessments(ts, n.noises, nassess)
    return Simulator(ts, n.model, scen, n.model.dynamics, getrealcost(n), realfinalcost)
end




################################################################################
# TODO: move Monte Carlo in another function
function simulate(simulator::Simulator, policy::AbstractPolicy)
    scenario = simulator.scenarios
    model = simulator.model

    # Get number of timesteps:
    ntime = size(scenario, 1)
    nsimu = size(scenario, 2)

    # Allocate
    stocks   = zeros(ntime, nsimu, model.dimStates)
    controls = zeros(ntime, nsimu, model.dimControls)
    costs    = zeros(nsimu)

    # Set first value of stocks equal to x0:
    for i in 1:nsimu
        stocks[1, i, :] = simulator.model.initialState
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

    return SimulationResult(costs, stocks, controls)
end
