################################################################################
# District.jl
################################################################################
# - Grid implements the global coordination problem.
# - It gathers nodes and interconnected edges
################################################################################


################################################################################
# Define whole microgrid with different interconnected noces.
abstract type AbstractGrid end

struct Grid <: AbstractGrid
    # Time span
    ts::AbstractTimeSpan
    # Nodes
    nodes::Vector{AbstractNode}
    # Edges
    net::AbstractNetwork
end
Grid(ts::AbstractTimeSpan) = Grid(ts, AbstractNode[], NoneNetwork())

nnodes(pb::Grid) = length(pb.nodes)
narcs(pb::Grid) = pb.net.narcs
ntimes(pb::Grid) = ntimesteps(pb.ts)


# Build models inside `grid` for decomposition
function build!(grid::Grid, xini::Dict, Interface::Type=PriceInterface;
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
function prodswap!(pb::Grid, mul::Vector{Float64})
    ntime = ntimes(pb) - 1
    for (id, d) in enumerate(pb.nodes)
        s1 = (id-1) * ntime + 1
        s2 = id * ntime
        swap!(d, mul[s1:s2])
    end
end
function flow!(pb::Grid, flow::Vector{Float64})
    ntime = ntimes(pb) - 1
    for (id, d) in enumerate(pb.nodes)
        s1 = (id-1) * ntime + 1
        s2 = id * ntime
        flow!(d, flow[s1:s2])
    end
end
# swap transport problem
transswap!(pb::Grid, mul::Vector{Float64}) = swap!(pb.net, mul)
# update graph exchange in nodes subproblems
function swap!(pb::Grid, mul::Vector{Float64})
    prodswap!(pb, mul)
    transswap!(pb, mul)
end

function swap!(pb::Grid, mul::Vector{Float64}, flow::Vector{Float64})
    # first, swap multipliers
    swap!(pb, mul)
    # then, swap current flows
    flow!(pb, flow)
end

# check consistency of grid with decomosition algorithm
function checkconsistency(pb::Grid, Interface::Type)
    chck1 = false ∉ isa.([m.conn for m in pb.nodes], Interface)
    return chck1
end

################################################################################
# Build global problem corresponding to Grid
################################################################################
initpos(pb::Grid) = vcat([h.model.initialState for h in pb.nodes]...)

function getproblem(pb::Grid, generation="reduction", nbins=10, noptscen=100)
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
    if generation == "reduction"
        laws = buildlaws(pb, noptscen, nbins)
    elseif generation == "total"
        # warning: usually intractable!!
        laws = Scenarios.prodprocess([towhitenoise(n.model.noises) for n in pb.nodes])
    end

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
function buildlink!(pb::Grid)
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


function builddynamic(pb::Grid)
    ntime = ntimes(pb)
    xindex = 1
    uindex = 1
    exdyn = Expr(:vect)

    params = Dict()
    params["text"] = loadweather(OutdoorTemperature(), pb.ts)

    for dev in pb.nodes
        # update irradiation in params
        pint, pext = get_irradiation(dev)
        params["pint"] = pint
        params["pext"] = pext

        # parse dynamics of nodes
        ex = parsebuilding(dev, xindex, uindex, dev.time.δt, params)
        xindex += nstocks(dev)
        uindex += ncontrols(dev)
        # push expression in global expression
        push!(exdyn.args, ex...)
    end

    return eval(:((t, x, u, w) -> $exdyn))
end

function buildcost(pb::Grid)
    function costgrid(m, t, x, u, w)
        # production cost
        cost = AffExpr(0.)
        xindex = 0
        uindex = 0
        for d in pb.nodes
            cost += objective(d, xindex, uindex)(m, t, x, u, w)
            xindex += nstocks(d)
            uindex += ncontrols(d)
        end
        # transport cost
        na = narcs(pb)
        # take abs |.| of flow q
        @variable(m, qp[1:na])
        @constraint(m, qp .>=  u[end-na+1:end])
        @constraint(m, qp .>= -u[end-na+1:end])
        # add transportation cost to cost
        cost += pb.net.k1*sum(qp) + pb.net.k2 * dot(qp, qp)
        return cost
    end
    return costgrid
end

# get index of flow importations in node
# remember: flow is last control in each node (u[end]) so its position
# is ncontrols
getflowindex(pb::Grid) = cumsum(ncontrols.(pb.nodes))

# Build coupling constraint Aq + F for a grid.
function buildconstr(pb::Grid)
    na = narcs(pb)
    A = pb.net.A
    findex = getflowindex(pb)

    constr(t, x, u, w) = A * u[end-na+1:end] + u[findex]
end


# Build probability laws for grid `pb`.
# WARNING
# Subject to curse of dimensionality (build laws in high dimension).
# TODO: can build very large 3D arrays...
function buildlaws(pb::Grid, nscen=100, nbins=10)
    # get total number of uncertainties
    nw = sum(nnoises.(pb.nodes))

    scenarios = zeros(Float64, ntimes(pb), nscen, nw)

    iw = 1

    for node in pb.nodes
        for ξ in node.noises
            ntw = nnoise(ξ)
            scenarios[:, :, iw:iw+ntw-1] = optscenarios(ξ, pb.ts, nscen)
            iw += ntw
        end
    end
    # then, we quantize vectors in `nw` dimensions
    return WhiteNoise(scenarios, nbins, KMeans())
end


# Build real cost
function getrealcost(pb::Grid)
    function realcost(t, x, u, w)
        cost = 0.
        for d in pb.nodes
            cost += getrealcost(d)(t, x, u, w)
        end
        return cost
    end
    return realcost
end

# TODO: clean definition of final cost (again!!!)
function buildfcost(pb::Grid)
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

function getrealfinalcost(pb::Grid)
    function fcost(x)
        cost = 0.
        xindex = 0
        for (id, d) in enumerate(pb.nodes)
            # get tank position
            postank = getposition(d, ElecHotWaterTank)
            # TODO: hardcoded!!!!
            cost += PENAL_TANK * (2. - x[postank+xindex])
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
