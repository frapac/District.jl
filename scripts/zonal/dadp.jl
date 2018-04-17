push!(LOAD_PATH, "..")

using District, StochDynamicProgramming
using Lbfgsb, Optim

include("problem.jl");
include("../plots/graph.jl")

srand(2713)

# Time span
ts = TimeSpan(180, 1)
# Noise discretization
nbins = 10
# Number of assessment in SDDP simulation
nassess = 1

# Build grid
pb, xini = twelvehouse(nbins=nbins)
nnodes = size(pb.nodes,1)

# Build SP problems in each node
build!(pb, xini, PriceInterface)

##################### SDDP #####
sim = Simulator(pb, nassess, generation="reduction", nbins=10)
params = District.get_sddp_solver()
params.max_iterations = 20
# Computing value functions 
sddp = solve_SDDP(sim.model, params, 2, 1)

# Simulating
pol = District.HereAndNowDP(sddp.bellmanfunctions)
resSDDP = District.simulate(sim, pol)
##################### SDDP #####

# Choosing weights for edges
q = getflow(pb,resSDDP.controls)
# Laplacian of incidence matrix
laplacian =  District.getlaplacian(pb.net.A, q)
# Get node assignments to clusters by spectral clustering
## /!\ Cluster number hardcoded
ncluster = 3 
membership = spectralclustering(laplacian, ncluster)

# Build zonal grid
pbreduced = District.reducegrid(pb, membership)

# Build SP problems in each zone
District.build!(pbreduced, xini, ZoneInterface, generation="reduction", nbins=10)


##################### DADP #####
algo = DADP(pbreduced, nsimu=1, nit=20)

# Initialization multiplier
mul0 = District.getinitialmultiplier(pbreduced)

# Compute optimal multipliers and value functions
f, grad! = District.oracle(pbreduced, algo)
gdsc = @time lbfgsb(f, grad!, mul0; iprint=1, pgtol=1e-5, factr=0., maxiter=20)

# Simulating
pol = District.DADPPolicy([algo.models[n.name].bellmanfunctions for n in pbreduced.nodes])
resDADP = District.simulate(sim, pol)
##################### DADP #####
