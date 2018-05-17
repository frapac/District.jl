function globallowerbound(pb::Grid, algo)

    lb = 0.
    for d in pb.nodes
        subpb = algo.models[d.name]
        lb += StochDynamicProgramming.lowerbound(subpb)
    end

    return lb + pb.net.cost
end

function elecload(res, sim, pb)
    # get expression of elecload
    ex = string(District.getelecload(pb))
    # parse expression and formate it
    ex = replace(ex, "[", "[:,:,")
    ex = replace(ex, "u", "res.controls")
    ex = replace(ex, "w", "sim.scenarios")
    # evaluation
    return eval(parse(ex))
end


hasmax(x) = x == :max
hasmax(x::Expr) = any(hasmax.(x.args))

function parsefcost(ex::Expr)
    hasmax(ex) && replacemax!(ex)
end

function replacemax!(ex::Expr)
    for pos in find(hasmax.(ex.args))
        subex = ex.args[pos]
        if subex.args[1] == :max
            ex.args[pos] = :z
        else
            replacemax!(ex.args[pos])
        end
    end
end

function addmax(ex)
    code = quote
        for i in 2:3
            @constraint($m, $z >= $(ex.args)[i])
        end
    end
    eval(code)
end

function recursivesampling(laws, sampler::DiscreteLawSampler, nbins, Δn)
    nwindows = ceil(Int, length(laws) / Δn)
    μtot = WhiteNoise[]

    ## we resample the noise sequentially
    i0, i1 = 1, Δn
    for nw in 1:nwindows
        μ = sampler(laws[i0:i1])
        push!(μtot, μ)
        i0 += Δn
        i1 = min(i1+Δn, length(laws))
    end

    return resample(μtot, nbins)
end
