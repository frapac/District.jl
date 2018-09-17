push!(LOAD_PATH, "..")

using District, StochDynamicProgramming
using Lbfgsb, Optim

include("../nodal/harsh.jl")
include("../plots/graph.jl")
include("../utils.jl")

# Initialize seed
srand(2713)
# Noise discretization
nbins = 10
# Build grid
grid, xini = twelvehouse(nbins=nbins)
N = District.nnodes(grid)
# Resampling of grid/zone model (for the SDDP model at each step of DADP)
resample = 5
# Maximum iteration for each algorithm
maxiterSDDP = 60
maxiterDADP = 20
# Number of Monte-Carlo simulation to estimate gradient (~ estimate import flows)
nmcsimu = 500


## Zonal decomposition
# Build SP problems in each node
build!(grid, xini, PriceInterface)

# Choosing weights for edges
q = vec(readcsv("results/zonal/$N/flows.csv"))

# Decompose grid and extract zone assignation
pb, membership = District.decomposegrid(grid, q, nclusters=4)

# DADP
# Rebuild SP problems and resampling noises from reference noises in each zone
build!(pb, xini, ZoneInterface, generation="reduction", nbins=resample)

# Initialization multiplier
mul0 = -District.getinitialmultiplier(pb)

algo = DADP(pb, nsimu=nmcsimu, nit=maxiterSDDP)

# Compute optimal multipliers and value functions
f, grad! = District.oracle(pb, algo)
gdsc = @time lbfgsb(f, grad!, mul0; iprint=1, pgtol=1e-5, factr=0., maxiter=maxiterDADP)
