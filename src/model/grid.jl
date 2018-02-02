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
function build!(grid::Grid, xini::Dict)
    for d in grid.nodes
        price = zeros(Float64, ntimes(grid)-1)
        # instantiate connection interface
        # TODO: dry connection size
        conn = PriceInterface(price, GraphConnection(6.))
        # add connection to particular device
        set!(d, conn)
        # build SP model inside `d`
        build!(d, xini[d])
    end
end

# update graph exchange in nodes subproblems
function swap!(pb::Grid, mul::Vector{Float64})
    ntime = ntimes(pb) - 1
    for (id, d) in enumerate(pb.nodes)
        s1 = (id-1) * ntime + 1
        s2 = id * ntime
        swap!(d, mul[s1:s2])
    end
    # swap transport problem
    swap!(pb.net, mul)
end


################################################################################
# Build global problem corresponding to Grid
################################################################################
initpos(pb::Grid) = vcat([h.model.initialState for h in pb.nodes]...)

function getproblem(pb::Grid)
    # to avoid world age problem, we rebuild elecload and gasload
    # right now.
    uindex = 1
    for d in pb.nodes
        # rebuild problem in Node `d`, starting from uindex
        buildload!(d, uindex)
        uindex += ncontrols(d)
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
    laws = buildlaws(pb)

    # TODO: how to define global fcost???
    spmodel = StochDynamicProgramming.LinearSPModel(ntimes(pb), ub,
                                                  x0, costm,
                                                  dynam,
                                                  tonoiselaws(laws),
                                                  info=:HD,
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
        # TODO: w index in builder
        ex = parsebuilding(dev, xindex, uindex, dev.time.δt, params)
        xindex += nstocks(dev)
        uindex += ncontrols(dev)
        # push expression in global expression
        push!(exdyn.args, ex...)
    end

    println(exdyn)
    return eval(:((t, x, u, w) -> $exdyn))
end

function buildcost(pb::Grid)
    function costgrid(m, t, x, u, w)
        cost = AffExpr(0.)
        for d in pb.nodes
            cost += objective(d)(m, t, x, u, w)
        end
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
