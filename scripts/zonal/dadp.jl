push!(LOAD_PATH, "..")

using District, StochDynamicProgramming
using Lbfgsb, Optim

include("problem.jl");
include("../plots/graph.jl")

srand(2713)

# Decomposition mode
mode = "ZONAL"
# Time span
ts = TimeSpan(180, 1)
# Noise discretization
nbins = 10
# Resampling of grid/zone model
resample = 10
# Number of assessment in SDDP simulation
nassess = 1
# Tolerance gap under which we consider the convergence of SDDP
tol = 0.1e-2
# Maximum iteration for each algorithm
maxiterSDDP = 30
maxiterDADP = 20
# Number of Monte-Carlo simulation to estimate gradient (~ estimate import flows)
nmcsimu = 1000

# Build grid
pb, xini = twelvehouse(nbins=nbins)
nnodes = size(pb.nodes,1)

# Build SP problems in each node
build!(pb, xini, PriceInterface)

if mode == "ZONAL"
    ##################### Zonal decomposition #####
    sim = Simulator(pb, nassess, generation="reduction", nbins=resample)
    params = District.get_sddp_solver()
    params.max_iterations = maxiterSDDP
    # Computing value functions 
    sddp = solve_SDDP(sim.model, params, 2, 1)

    # Simulating
    pol = District.HereAndNowDP(sddp.bellmanfunctions)
    resSDDP = District.simulate(sim, pol)
    costs, states, controls = StochDynamicProgramming.simulate(sddp, nassess)
    lb = sddp.stats.lower_bounds
    gap = costs./lb - 1
    conv = gap .< tol
    nit = min(maxiterSDDP, size(conv[.!conv],1) + 1)
    convergencetime = sum(sddp.stats.exectime[1:nit])

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

    # ##################### Zonal decomposition #####

    # # Rebuild SP problems and resampling noises from scratch in each zone
    # District.build!(pbreduced, xini, ZoneInterface, generation="reduction", nbins=resample)

    # # Initialization multiplier
    # mul0 = -District.getinitialmultiplier(pbreduced)

    # algo = DADP(pbreduced, nsimu=nmcsimu, nit=maxiterSDDP)

    # # Compute optimal multipliers and value functions
    # f, grad! = District.oracle(pbreduced, algo)

    # gdsc = @time lbfgsb(f, grad!, mul0; iprint=1, pgtol=1e-5, factr=0., maxiter=maxiterDADP)

    # Simulating
    #pol = District.DADPPolicy([algo.models[n.name].bellmanfunctions for n in pbreduced.nodes])
    #resDADP = District.simulate(sim, pol)

elseif mode == "NODAL"
    # Initialization multiplier
    mul0 = -District.getinitialmultiplier(pb)
    
    # Algorithm to resolve value functions problem
    algo = DADP(pb, nsimu=nsimu, nit=maxiterSDDP)

    # Compute optimal multipliers and value functions
    f, grad! = District.oracle(pb, algo)

    gdsc = @time lbfgsb(f, grad!, mul0; iprint=1, pgtol=1e-5, factr=0., maxiter=maxiterDADP)

    # Simulating
    #pol = District.DADPPolicy([algo.models[n.name].bellmanfunctions for n in pb.nodes])
    #resDADP = District.simulate(sim, pol)
end



