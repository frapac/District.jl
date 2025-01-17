################################################################################
# District.jl
################################################################################
# Implement Monte-Carlo simulation.
# - Simulator stores all information needed for simulation.
# - Simulation returns a SimulationResult object.
################################################################################

export Simulator, SimulationResult


################################################################################
# Simulation results
"""
    struct SimulationResult
        # costs along assessment scenarios. size = (nscen,)
        costs::Vector{Float64}
        # stocks along assessment scenarios. size = (ntime, nscen, nx)
        stocks::Array{Float64, 3}
        # controls along assessment scenarios. size = (ntime, nscen, nu)
        controls::Array{Float64, 3}
    end

Object to store results of a simulation.
"""
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
"""
    struct Simulator
        # names of states, controls and uncertainties
        names::Dict
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

Simulator object.


# Construction

    Simulator(n::AbstractNode, nassess::Int)

Build Simulator corresponding to Node `n` and to `nassess` scenarios.

    Simulator(pb::AbstractGrid, nassess::Int)

Build Simulator corresponding to Grid `pb` and to `nassess` scenarios.
"""
struct Simulator
    # names of states, controls and uncertainties
    names::Dict
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
function Simulator(n::AbstractNode, nassess::Int, gen::AbstractScenarioGenerator=OutSampleScenarios())
    xname, uname, wname = getcorrespondance(n)
    names = Dict(:x=>xname, :u=>uname, :w=>wname)
    ts = n.time
    scen = generate(gen, n, nassess)
    return Simulator(names, ts, n.model, scen, n.model.dynamics, getrealcost(n), realfinalcost)
end


# adapt Simulator for grid
function Simulator(pb::AbstractGrid, nassess::Int;
                   sampler::DiscreteLawSampler=DiscreteLawSampler(1, 1, 1),
                   gen::AbstractScenarioGenerator=OutSampleScenarios())
    # Get names label
    xname, uname, wname = getcorrespondance(pb)
    names = Dict(:x=>xname, :u=>uname, :w=>wname)
    # Generate time span
    ts = pb.ts

    # Model generation
    # if specified, set random seed
    (gen.seed > 0) && srand(gen.seed)
    # build global problem
    spmodel = getproblem(pb, sampler)
    # generate scenarios
    scen = generate(gen, pb, nassess)
    return Simulator(names, ts, spmodel, scen, spmodel.dynamics,
                     getrealcost(pb), getrealfinalcost(pb))
end

function show(io::IO, sim::Simulator)
    println("* Simulation period:  day $(sim.ts.day) (ndays: $(sim.ts.ndays))")
    println("* Number of scenarios: ", size(sim.scenarios, 2))
end

################################################################################
# TODO: move Monte Carlo in another function
"""
    simulate(sim::Simulator, policy::AbstractPolicy)

Simulate strategies specified by `policy` with `sim` Simulator.
Return a `SimulationResult` object.
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
            # TODO
            costs[k] += simulator.realcost(t, x, u, ξ)

            m = policy.problem
            #Sum over instantaneous costs only
            #= costs[k] += getobjectivevalue(m) - getvalue(m[:alpha]) =#
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


################################################################################
# UTILS
################################################################################
function getcorrespondance(h::AbstractNode)
    xname = String[]
    uname = String[]
    wname = String[]

    for n in h.devices
        nx = nstates(n)
        nu = ncontrols(n)
        push!(xname, fill(getname(n), nx)...)
        push!(uname, fill(getname(n), nu)...)
    end

    for w in h.noises
        nw = nnoise(w)
        push!(wname, fill(getname(w), nw)...)
    end

    return xname, uname, wname
end

function getcorrespondance(pb::AbstractNodalGrid; lastzoneindex=0)
    xname = String[]
    uname = String[]
    wname = String[]

    for (iin, n) in enumerate(pb.nodes)
        xn, un, wn = getcorrespondance(n)
        push!(xname, ("Node $(lastzoneindex+iin): " .* xn)...)
        push!(uname, ("Node $(lastzoneindex+iin): " .* un)...)
        push!(wname, ("Node $(lastzoneindex+iin): " .* wn)...)
    end

    for edge in 1:narcs(pb)
        pos = find(x->(x!=0), pb.net.A[:, edge])
        # one edge joins only two nodes
        @assert length(pos) == 2
        push!(uname, ("Connection $(lastzoneindex+pos[1]) <--> $(lastzoneindex+pos[2])"))
    end

    ninj = ninjection(pb)
    if ninj > 0
        for inj in 1:ninj
            push!(uname, ("Ext. Import $inj" ))
        end
    end


    return xname, uname, wname
end

function getcorrespondance(pb::ZonalGrid)
    xname = String[]
    uname = String[]
    wname = String[]

    lastzoneindex = 0
    for (iin, n) in enumerate(pb.nodes)
        xn, un, wn = getcorrespondance(n, lastzoneindex=lastzoneindex)
        xname = vcat(xname, xn)
        uname = vcat(uname, un)
        wname = vcat(wname, wn)
        lastzoneindex+=length(pb.nodes)
    end

    return xname, uname, wname
end


"""
    getlabel(sim::Simulator, k::Symbol)

Print labels of objects inside Simulator `sim`.

Different choices of `k` are:
* `:x`: print states labels ;
* `:u`: print controls labels ;
* `:w`: print noises labels .

"""
function getlabel(sim::Simulator, k::Symbol)
    names = sim.names[k]
    println("="^30)
    for (iix, n) in enumerate(names)
        println("$k[$iix]: $n")
    end
    println("="^30)
end

function getlabel(pb::AbstractGrid, k::Symbol)
    xname, uname, wname = getcorrespondance(pb)
    namesdict = Dict(:x=>xname, :u=>uname, :w=>wname)
    names = namesdict[k]
    println("="^30)
    for (iix, n) in enumerate(names)
        println("$k[$iix]: $n")
    end
    println("="^30)
end
