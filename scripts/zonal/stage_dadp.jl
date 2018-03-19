push!(LOAD_PATH, "..")

using District, StochDynamicProgramming
using Lbfgsb, Optim

srand(2713)

ALGO = "SDDP"


# Time span
ts = TimeSpan(180, 1)

# Noise discretization
nbins = 10

# Build graph
#include("toygraph.jl")
pb, xini = twelvehouse(nbins=nbins)

nnodes = size(pb.nodes,1)

# Build SP problems in each node
build!(pb, xini, PriceInterface)



# Launch Gradient Descent!
#= x0 = zeros(Float64, size(A, 1)*95) =#
p = EDFPrice(ts).price[1:end-1]
# initialization multiplicator
mul0 = p
for i in 1:nnodes-1
    mul0 = vcat(mul0, p)
end

##################### SDDP #####
    nassess = 1
    sim = Simulator(pb, nassess, generation="reduction", nbins=30)
    params = District.get_sddp_solver()
    params.max_iterations = 50
    sddp = solve_SDDP(sim.model, params, 2, 1)

    pol = District.HereAndNowDP(sddp.bellmanfunctions)
    res = District.simulate(sim, pol)
##################### SDDP #####

q = mean(mean(abs.(getflow(res.controls)),2),1)

laplacian =  getlaplacian(pb.net.A , q[1,1,:])

res = spectralclustering(laplacian, ncluster)


# Coloring nodes
membership = res.assignments

for h in 1:nnodes
    if membership[h] == 1
        vec1 = hcat(vec1, pb.nodes[h])
    elseif membership[h] == 2
        vec2 = hcat(vec2, pb.nodes[h])
    elseif membership[h] == 3
        vec3 = hcat(vec3, pb.nodes[h])
    end
end



# Build zones 
z1 = Zone(ts, vec1, net1)
z2 = Zone(ts, vec2)
z3 = Zone(ts, vec3)

pbreduced = Grid(ts, [z1, z2, z3], netreduced)

# Build SP problems in each zone
zonebuild!(pbreduced, xini, PriceInterface)


# Select algorithm to solve global problem
algo = DADP(pbreduced, nsimu=1, nit=10)

f, grad! = District.oracle(pbreduced, algo)
gdsc = @time lbfgsb(f, grad!, mul0; iprint=1, pgtol=1e-5, factr=0., maxiter=40)

pol = District.DADPPolicy([algo.models[n.name].bellmanfunctions for n in pbreduced.nodes])
res = District.simulate(sim, pol, )

