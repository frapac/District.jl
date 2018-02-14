################################################################################
# District.jl
################################################################################
# Implement Monte-Carlo simulation.
# - Simulator stores all information needed for simulation.
# - Simulation returns a SimulationResult object.
################################################################################

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
# overload Base.show
function show(io::IO, res::SimulationResult)
    @printf("Costs: %.4f ± %.4f €\n", mean(.25res.costs), std(.25res.costs))
end

################################################################################
# Simulator
struct Simulator
    # time period considered
    ts::TimeSpan
    # SP Model to simulate
    model::StochDynamicProgramming.SPModel
    # assessment scenarios
    scenarios::Array{Float64, 3}
    # simulation's dynamics
    realdynamic::Function
    # simulations's costs
    realcost::Function
    # simulation's final cost
    realfinalcost::Function
end


# TODO: currently dynamics is the same as in optimization model
"""
    Simulator(n::AbstractNode, nassess::Int)

Build Simulator corresponding to Node `n` and to `nassess` scenarios.
"""
function Simulator(n::AbstractNode, nassess::Int)
    ts = n.time
    scen = genassessments(ts, n.noises, nassess)
    return Simulator(ts, n.model, scen, n.model.dynamics, getrealcost(n), realfinalcost)
end


# adapt Simulator for grid
function Simulator(pb::AbstractGrid, nassess::Int;
                   generation="reduction", nbins=10, noptscen=100)
    ts = pb.ts
    scen = genassessments(pb, nassess)
    # build global problem
    spmodel = getproblem(pb, generation, nbins, noptscen)
    return Simulator(ts, spmodel, scen, spmodel.dynamics,
                     getrealcost(pb), getrealfinalcost(pb))
end



################################################################################
# TODO: move Monte Carlo in another function
"""
    simulate(sim::Simulator, policy::AbstractPolicy)

Simulate strategies specified by `policy` on `sim` Simulator.
Return a SimulationResult object.
"""
function simulate(simulator::Simulator, policy::AbstractPolicy)
    scenario = simulator.scenarios
    model = simulator.model

    # Get number of timesteps and number of simulations
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

    # simulate
    @showprogress for t=1:ntime-1
        # update problem inside policy
        buildproblem!(policy, simulator.model, t)

        for k in 1:nsimu
            # get previous state:
            x = stocks[t, k, :]
            # and current noise
            ξ = scenario[t, k, :]

            # compute decisions with policy:
            u = policy(x, ξ)
            # get future state with real dynamics
            xf = simulator.realdynamic(t, x, u, ξ)

            # update
            costs[k] += simulator.realcost(t, x, u, ξ)
            stocks[t+1, k, :] = xf
            controls[t, k, :] = u
        end
    end
    # compute final costs along scenarios
    for k = 1:nsimu
        costs[k] += simulator.realfinalcost(stocks[end, k, :])
    end

    return SimulationResult(costs, stocks, controls)
end
