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
res = District.simulate(sim, pol)
##################### SDDP #####


# Time-scenario mean flow on edges
q = getflow(res.controls)

# Laplacian of incidence matrix
laplacian =  getlaplacian(pb.net.A , q)

# Apply spectral clustering algorithm
## Cluster number hardcoded
ncluster = 3 
clusterresult = spectralclustering(laplacian, ncluster)

# Get node assignments to clusters
membership = clusterresult.assignments 

# Fill zones
zones, netreduced = District.fillzones(pb, membership)

# Build zonal grid
pbreduced = District.ZonalGrid(ts, zones, netreduced)

xname, uname, wname = District.getcorrespondance(pbreduced)
names = Dict(:x=>xname, :u=>uname, :w=>wname)

# Build SP problems in each zone
District.zonebuild!(pbreduced, xini, PriceInterface, generation="reduction", nbins=1)


##################### DADP #####
algo = DADP(pbreduced, nsimu=1, nit=5)


p = EDFPrice(ts).price[1:end-1]

# Initialization multiplicator
mul0 = p
nbordernodes = sum(length.([zone.bordernodes for zone in zones]))
for i in 1:nbordernodes-1
    mul0 = vcat(mul0, p)
end

# Compute f:mul and grad!:(mul, storage)
f, grad! = District.oracle(pbreduced, algo)
gdsc = @time lbfgsb(f, grad!, mul0; iprint=1, pgtol=1e-5, factr=0., maxiter=10)

pol = District.DADPPolicy([algo.models[n.name].bellmanfunctions for n in pbreduced.nodes])
#res = District.simulate(sim, pol, )
##################### DADP #####
