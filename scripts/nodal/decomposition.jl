################################################################################
# District.jl
################################################################################
# Generate academic problem and solve SP model by decomposition.
################################################################################

push!(LOAD_PATH, "..")


using District
using StochDynamicProgramming
using Lbfgsb, Ipopt
include("harsh.jl")
include("solvers.jl")

srand(2713)

ALGO = "DADP"

pb, xini = sixhouse(nbins=10)

if ALGO == "DADP"
    build!(pb, xini, PriceInterface, maxflow=6.)
    _, algo = bfgs(pb, nsimu=1000)
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
    _, algo = quantdec(pb, nsimu=500)
    pol = District.DADPPolicy([algo.models[n.name].bellmanfunctions for n in pb.nodes])
elseif ALGO == "SDDP"
    build!(pb, xini, PriceInterface, maxflow=6.)
    sim  = Simulator(pb, 1, generation="reduction", nbins=80, outsample=false)
    algo = runsddp(sim.model)
end
