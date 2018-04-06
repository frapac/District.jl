push!(LOAD_PATH, "..")

using District, StochDynamicProgramming
using Lbfgsb, Optim

srand(2713)

include("../graph.jl")
include("../problem.jl")

# Time span
ts = TimeSpan(180, 1)

# Noise discretization
nbins = 1

# Build grid
pb, xini = twelvehouse(nbins=nbins)

nnodes = size(pb.nodes,1)

# Build SP problems in each node
build!(pb, xini, PriceInterface)

##################### SDDP #####
nassess = 1
sim = Simulator(pb, nassess, generation="reduction", nbins=1)
params = District.get_sddp_solver()
params.max_iterations = 5
# Computing value functions 
sddp = solve_SDDP(sim.model, params, 2, 1)

pol = District.HereAndNowDP(sddp.bellmanfunctions)
resSDDP = District.simulate(sim, pol)
##################### SDDP #####


# Choice of weights for edges
q = getflow(resSDDP.controls)

# Get node assignments to clusters by spectral clustering
## Cluster number hardcoded
ncluster = 3 
membership = spectralclustering(pb.net.A, q, ncluster)

# membership = vec([1 1 1 2 2 2 3 3 3 3 3 3])

# Build zonal grid
pbreduced = District.reducegrid(pb, membership)

# Build SP problems in each zone
District.build!(pbreduced, xini, ZoneInterface, generation="reduction", nbins=1)


##################### DADP #####
algo = DADP(pbreduced, nsimu=1, nit=5)

# Initialization multiplier
p = EDFPrice(ts).price[1:end-1]
sizelambda = sum(District.nbordernodes.([zone for zone in pbreduced.nodes]))
mul0 = District.getinitialmultiplier(p, sizelambda)

# Compute optimal multipliers
f, grad! = District.oracle(pbreduced, algo)
gdsc = @time lbfgsb(f, grad!, mul0; iprint=1, pgtol=1e-5, factr=0., maxiter=10)

# Once optimal multipliers obtained, define optimal policy
pol = District.DADPPolicy([algo.models[n.name].bellmanfunctions for n in pbreduced.nodes])
resDADP = District.simulate(sim, pol)
##################### DADP #####
