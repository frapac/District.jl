################################################################################
# District.jl
################################################################################
# - Grid implements the global coordination problem.
# - It gathers nodes and interconnected edges
################################################################################


################################################################################
# Define whole microgrid with different interconnected noces.
abstract type AbstractGrid end

abstract type AbstractNodalGrid <: AbstractGrid end

struct Grid <: AbstractNodalGrid
    # Time span
    ts::AbstractTimeSpan
    # Nodes
    nodes::Vector{AbstractNode}
    # Edges
    net::AbstractNetwork
end
Grid(ts::AbstractTimeSpan) = Grid(ts, AbstractNode[], NoneNetwork())



nnodes(pb::AbstractNodalGrid) = length(pb.nodes)
narcs(pb::AbstractGrid) = pb.net.narcs
ntimes(pb::AbstractGrid) = ntimesteps(pb.ts)
ninjection(pb::Grid) = 0


# Build models inside `grid` for decomposition
function build!(grid::AbstractNodalGrid, xini::Dict, Interface::Type=PriceInterface;
                maxflow=6., tau=1.)
    for d in grid.nodes
        price = zeros(Float64, ntimes(grid)-1)
        # instantiate connection interface
        if Interface == QuadInterface
            flow = zeros(Float64, ntimes(grid)-1)
            conn = Interface(tau, price, flow, GraphConnection(maxflow))
        else
            conn = Interface(price, GraphConnection(maxflow))
        end

        # add connection to particular device
        set!(d, conn)
        # build SP model inside `d`
        build!(d, xini[d])
    end
end

# swap production problem
function prodswap!(pb::AbstractGrid, mul::Vector{Float64})
    slambda = 0
    for d in pb.nodes
        s1 = slambda + 1
        slambda += connectionsize(d)
        s2 = slambda
        swap!(d, mul[s1:s2])
    end
end

function flow!(pb::AbstractGrid, flow::Vector{Float64})
    ntime = ntimes(pb) - 1
    for (id, d) in enumerate(pb.nodes)
        s1 = (id-1) * ntime + 1
        s2 = id * ntime
        flow!(d, flow[s1:s2])
    end
end

# swap transport problem
transswap!(pb::AbstractGrid, mul::Vector{Float64}) = swap!(pb.net, mul)
# update graph exchange in nodes subproblems
function swap!(pb::AbstractGrid, mul::Vector{Float64})
    prodswap!(pb, mul)
    transswap!(pb, mul)
end

function swap!(pb::AbstractGrid, mul::Vector{Float64}, flow::Vector{Float64})
    # first, swap multipliers
    swap!(pb, mul)
    # then, swap current flows
    flow!(pb, flow)
end

# check consistency of grid with decomosition algorithm
function checkconsistency(pb::AbstractGrid, Interface::Type)
    chck1 = false ∉ isa.([m.conn for m in pb.nodes], Interface)
    return chck1
end

################################################################################
# Build global problem corresponding to Grid
################################################################################
initpos(pb::AbstractNodalGrid) = vcat([h.model.initialState for h in pb.nodes]...)

function getproblem(pb::AbstractNodalGrid, sampler::DiscreteLawSampler)
    # to avoid world age problem, we rebuild elecload and gasload
    # right now.
    uindex = 1
    windex = 1
    for d in pb.nodes
        # rebuild problem in Node `d`, starting from uindex
        buildload!(d, uindex, windex)
        uindex += ncontrols(d)
        windex += nnoises(d)
    end


    # get global initial position
    x0 = initpos(pb)
    xb = vcat(xbounds.(pb.nodes)...)
    ub = vcat(ubounds.(pb.nodes)...)
    # add flow control
    ub = vcat(ub, [(-m, m) for m in pb.net.maxflow])

    # TO DO : choose proper bounds
    for i in 1:ninjection(pb)
        ub = vcat(ub, [pb.bounds[i]])
    end

    # In order ...
    # reset and rebuild linkage in nodes (mainly update Expr junction)
    buildlink!(pb)
    # build global dynamics
    dynam = builddynamic(pb)
    # build global cost
    costm = buildcost(pb)
    # build coupling constraint Aq+f, modeling graph topology
    constr = buildconstr(pb)
    # build global noise laws (warning: |WW| may be very large)
    laws = sampler([towhitenoise(n.model.noises) for n in pb.nodes])

    spmodel = StochDynamicProgramming.LinearSPModel(ntimes(pb), ub,
                                                  x0, costm,
                                                  dynam,
                                                  tonoiselaws(laws),
                                                  info=:HD,
                                                  Vfinal=buildfcost(pb),
                                                  eqconstr=constr)

    set_state_bounds(spmodel, xb)
    return spmodel
end

# TODO: possible side effect with reset!
# Only call this function at the end!!!
function buildlink!(pb::AbstractNodalGrid)
    uindex = 0
    windex = 0
    for node in pb.nodes
        # reset all devices Expr
        reset!.(node.devices)
        # build link with global index
        buildlink!(node, uindex, windex)
        # update global index
        uindex += ncontrols(node)
        windex += nnoises(node)
    end
end


function builddynamic(pb::AbstractNodalGrid)
    ntime = ntimes(pb)
    xindex = 1
    uindex = 1
    exdyn = Expr(:vect)

    params = Dict()
    params["text"] = loadweather(OutdoorTemperature(), pb.ts)

    for dev in pb.nodes
        # update irradiation in params
        if hasdevice(dev, R6C2)
            pint, pext = get_irradiation(dev)
            params["pint"] = pint
            params["pext"] = pext
        end

        # parse dynamics of nodes
        ex = parsebuilding(dev, xindex, uindex, dev.time.δt, params)
        xindex += nstocks(dev)
        uindex += ncontrols(dev)
        # push expression in global expression
        push!(exdyn.args, ex...)
    end

    return eval(:((t, x, u, w) -> $exdyn))
end

function buildcost(pb::AbstractNodalGrid)
    ninj = ninjection(pb)
    function costgrid(m, t, x, u, w)
        # production cost
        cost = AffExpr(0.)
        xindex = 0
        uindex = 0
        # Instantaneous cost of all nodes
        for d in pb.nodes
            cost += objective(d, xindex, uindex)(m, t, x, u, w)
            xindex += nstocks(d)
            uindex += ncontrols(d)
        end
        # transport cost
        na = narcs(pb)
        # take abs |.| of flow q
        @variable(m, qp[1:na])
        @constraint(m, qp .>=  u[end-ninj-na+1:end-ninj])
        @constraint(m, qp .>= -u[end-ninj-na+1:end-ninj])
        # add transportation cost to cost
        cost += pb.net.k1*sum(qp) + pb.net.k2 * dot(qp, qp)

        # Zone connection cost
        if ninj > 0
            cost += objective(pb)(m, t, x, u, w)
        end
        return cost
    end
    return costgrid
end

# get index of flow importations in node
# remember: flow is last control in each node (u[end]) so its position
# is ncontrols
getflowindex(pb::AbstractNodalGrid) = cumsum(ncontrols.(pb.nodes))


# Build coupling constraint Aq + F for a grid.
function buildconstr(pb::AbstractNodalGrid)
    na = narcs(pb)
    ninj = ninjection(pb)
    A = pb.net.A
    findex = getflowindex(pb)

    function constr(t, x, u, w)
        cstr = A * u[end-ninj-na+1:end-ninj] + u[findex]
        if ninj > 0
            cstr += pb.conn.linker * u[end-ninj+1:end]
        end
        return cstr
    end
    return constr
end

# Build real cost
function getrealcost(pb::AbstractNodalGrid)
    warn("Building real cost for `Grid` is currently broken")
    function realcost(t, x, u, w)
        cost = 0.
        for d in pb.nodes
            cost += getrealcost(d)(t, x, u, w)
        end
        # add network cost
        na = narcs(pb)
        cost += getcost(pb.net, u[end-na+1:end])
        return cost
    end
    return realcost
end

# TODO: clean definition of final cost (again!!!)
function buildfcost(pb::AbstractNodalGrid)
    function final_cost(model, m)
        alpha = m[:alpha]
        x = m[:x]
        u = m[:u]
        xf = m[:xf]

        xindex = 0
        @JuMP.variable(m, z1[1:nnodes(pb)]>=0)
        for (id, d) in enumerate(pb.nodes)
            # get tank position
            postank = getposition(d, ElecHotWaterTank)
            @JuMP.constraint(m, z1[id] >= PENAL_TANK * (2. - xf[postank+xindex]))
            xindex += nstocks(d)
        end

        @JuMP.constraint(m, alpha == sum(z1))
    end
    return final_cost
end

# return elec load of grid as expression
function getelecload(pb::Grid)
    ex = Expr(:call, :+)
    uindex, windex = (1, 1)
    for (id, d) in enumerate(pb.nodes)
        push!(ex.args, getload(d, uindex, windex))
        uindex += ncontrols(d)
        windex += nnoises(d)
    end
    return ex
end

function getrealfinalcost(pb::AbstractNodalGrid)
    function fcost(x)
        cost = 0.
        xindex = 0
        for (id, d) in enumerate(pb.nodes)
            # get tank position
            postank = getposition(d, ElecHotWaterTank)
            # TODO: hardcoded!!!!
            cost += PENAL_TANK * max(2. - x[postank+xindex], 0.)
            xindex += nstocks(d)
        end
        return cost
    end
    return fcost
end

################################################################################
# PRINT
################################################################################
function show(io::IO, pb::Grid)
    println("Grid with ", nnodes(pb), " nodes.")
    for n in pb.nodes
        println("- Node ", n.name, " (n stocks: ", nstocks(n),")")
    end
end
