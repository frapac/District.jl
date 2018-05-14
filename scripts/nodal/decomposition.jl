################################################################################
# District.jl
################################################################################
# Generate academic problem and solve SP model by decomposition.
################################################################################

push!(LOAD_PATH, "..")


using District, Scenarios
using StochDynamicProgramming
using Lbfgsb, Ipopt
include("harsh.jl")
include("solvers.jl")
include("../utils.jl")

srand(2713)

ALGO = "DADP"

pb, xini = house48(nbins=10)

if ALGO == "DADP"
    build!(pb, xini, PriceInterface, maxflow=6.)
    _, algo = bfgs(pb, nsimu=500, sddpit=20, nit=50)
    pol = District.DADPPolicy([algo.models[n.name].bellmanfunctions for n in pb.nodes])
elseif ALGO == "QADP"
    build!(pb, xini, FlowInterface, maxflow=6.)
    _, algo = qadp(pb, nsimu=1000)
    pol = District.DADPPolicy([algo.models[n.name].bellmanfunctions for n in pb.nodes])
elseif ALGO == "IPOPT"
    build!(pb, xini, PriceInterface, maxflow=6.)
    _, algo = ipopt(pb)
    pol = District.DADPPolicy([algo.models[n.name].bellmanfunctions for n in pb.nodes])
elseif ALGO == "PADP"
    build!(pb, xini, FlowInterface, maxflow=6.)
    _, algo = quantdec(pb, nsimu=1000)
    pol = District.DADPPolicy([algo.models[n.name].bellmanfunctions for n in pb.nodes])
elseif ALGO == "SDDP"
    build!(pb, xini, PriceInterface, maxflow=6.)
    model = @time District.getproblem(pb, DiscreteLawSampler(1, 1, Δn=1))
    laws = [District.towhitenoise(n.model.noises) for n in pb.nodes]
    μquant = recursivesampling(laws, DiscreteLawSampler(10, 5, Δn=4), 100, 12)
    model.noises = District.tonoiselaws(μquant)
    algo = runsddp(model, nit=1500, ncuts=200)
end
