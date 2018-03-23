push!(LOAD_PATH, "..")

using District, StochDynamicProgramming
using Lbfgsb, Optim

srand(2713)

include("../graph.jl")

# Time span
ts = TimeSpan(180, 1)

# Noise discretization
nbins = 10

# Build grid
#include("toygraph.jl")
include("../problem.jl")
pb, xini = twelvehouse(nbins=nbins)

nnodes = size(pb.nodes,1)

# Build SP problems in each node
build!(pb, xini, PriceInterface)




##################### SDDP #####
    nassess = 1
    sim = Simulator(pb, nassess, generation="reduction", nbins=30)
    params = District.get_sddp_solver()
    params.max_iterations = 1
    sddp = solve_SDDP(sim.model, params, 2, 1)

    pol = District.HereAndNowDP(sddp.bellmanfunctions)
    res = District.simulate(sim, pol)
##################### SDDP #####




# Time-scenario mean flow on edges
q = mean(mean(abs.(getflow(res.controls)),2),1)

# Laplacian of incidence matrix
laplacian =  getlaplacian(pb.net.A , q[1,1,:])

# Apply spectral clustering algorithm
## Cluster number hardcoded
ncluster = 3 
clusterresult = spectralclustering(laplacian, ncluster)

# Get node assignments to clusters
membership = clusterresult.assignments 

# Fill zones
zones, netreduced = fillzones(pb, membership)

# Build zonal grid
pbreduced = ZonalGrid(ts, zones, netreduced)

# Build SP problems in each zone
zonebuild!(pbreduced, xini, PriceInterface)




##################### DADP #####
	algo = DADP(pbreduced, nsimu=1, nit=10)

	p = EDFPrice(ts).price[1:end-1]

	# Initialization multiplicator
	mul0 = p
	for i in 1:nnodes-1
	    mul0 = vcat(mul0, p)
	end

	# Compute f:mul and grad!:(mul, storage)
	f, grad! = District.oracle(pbreduced, algo)
	gdsc = @time lbfgsb(f, grad!, mul0; iprint=1, pgtol=1e-5, factr=0., maxiter=40)

	pol = District.DADPPolicy([algo.models[n.name].bellmanfunctions for n in pbreduced.nodes])
	res = District.simulate(sim, pol, )
##################### DADP #####
