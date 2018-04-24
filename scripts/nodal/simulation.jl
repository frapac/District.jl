################################################################################
# District.jl
################################################################################
# Scripts to simulate policies with previously computed cuts.
################################################################################

push!(LOAD_PATH, "..")


using District, Scenarios
using StochDynamicProgramming
include("harsh.jl")

# Choose algorithm
ALGO = "PADP"

# fix seed for reproductability
srand(2713)

# Build problem
pb, xini = threehouse(nbins=10)
# Get number of nodes
N = District.nnodes(pb)
# Build SP model inside nodes
build!(pb, xini, PriceInterface, maxflow=6.)
# Build corresponding simulator
sim = Simulator(pb, 5000,
                sampler=DiscreteLawSampler(1, 1, 1),
                gen=InSampleScenarios())

FOLD = "$N-d"
COM = ".1000"
# Generate policies
if ALGO == "DADP"
    vfs = [StochDynamicProgramming.read_polyhedral_functions("results/nodal/$FOLD/dadp$COM/$i-cuts.csv") for i in 1:N]
    pol = District.DADPPolicy(vfs)
elseif ALGO == "PADP"
    vfs = [StochDynamicProgramming.read_polyhedral_functions("results/nodal/$FOLD/qadp$COM/$i-cuts.csv") for i in 1:N]
    pol = District.DADPPolicy(vfs)
elseif ALGO == "SDDP"
    vfs = StochDynamicProgramming.read_polyhedral_functions("results/nodal/$FOLD/sddp$COM/cuts.csv")
    pol = District.HereAndNowDP(vfs)
end

# Simulate!
res = District.simulate(sim, pol)
