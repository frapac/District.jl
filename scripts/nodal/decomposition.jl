################################################################################
# District.jl
################################################################################
# Generate academic problem and solve SP model by decomposition.
################################################################################

push!(LOAD_PATH, "..")


using District, Scenarios
using StochDynamicProgramming
using Lbfgsb #, Ipopt
include("harsh.jl")
include("solvers.jl")

srand(2713)

ALGO = "DADP"

pb, xini = threehouse(nbins=10)

if ALGO == "DADP"
    build!(pb, xini, PriceInterface, maxflow=6.)
    _, algo = bfgs(pb, nsimu=500, sddpit=40, nit=20)
elseif ALGO == "QADP"
    build!(pb, xini, FlowInterface, maxflow=6.)
    _, algo = qadp(pb, nsimu=1000)
elseif ALGO == "IPOPT"
    build!(pb, xini, PriceInterface, maxflow=6.)
    prob, algo, vals, histmul = ipopt(pb, nsimu=500, nit=25)
elseif ALGO == "PADP"
    build!(pb, xini, FlowInterface, maxflow=6.)
    prob, algo, vals, histmul = quantdec(pb, nsimu=500)
elseif ALGO == "SDDP"
    build!(pb, xini, PriceInterface, maxflow=6.)
    model = @time District.getproblem(pb, DiscreteLawSampler(10, 5, 1))
    algo = runsddp(model)
end
