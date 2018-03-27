################################################################################
# District.jl
################################################################################
# Scripts to simulate policies with previously computed cuts.
################################################################################

push!(LOAD_PATH, "..")


using District
using StochDynamicProgramming
include("problem.jl")

# Choose algorithm
ALGO = "SDDP"

# fix seed for reproductability
srand(2713)

# Build problem
pb, xini = twelvehouse(nbins=10)
# Get number of nodes
N = District.nnodes(pb)
# Build SP model inside nodes
build!(pb, xini, PriceInterface, maxflow=6.)
# Build corresponding simulator
sim = Simulator(pb, 100, generation="reduction", nbins=50, outsample=true)

# Generate policies
if ALGO == "DADP"
    vfs = [StochDynamicProgramming.read_polyhedral_functions("results/nodal/$N/dadp/$i-cuts.csv") for i in 1:N]
    pol = District.DADPPolicy(vfs)
elseif ALGO == "PADP"
    vfs = [StochDynamicProgramming.read_polyhedral_functions("results/nodal/$N/qadp/$i-cuts.csv") for i in 1:N]
    pol = District.DADPPolicy(vfs)
elseif ALGO == "SDDP"
    vfs = StochDynamicProgramming.read_polyhedral_functions("results/nodal/$N/sddp/cuts.csv")
    pol = District.HereAndNowDP(vfs)
end

# Simulate!
res = District.simulate(sim, pol)
