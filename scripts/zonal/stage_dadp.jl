
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
# Select algorithm to solve global problem
algo = DADP(pb, nsimu=1, nit=10)


# Launch Gradient Descent!
#= x0 = zeros(Float64, size(A, 1)*95) =#
p = EDFPrice(ts).price[1:end-1]
# initialization multiplicator
x0 = p
for i in 1:nnodes-1
    x0 = vcat(x0, p)
end

if ALGO == "SDDP"
    # TODO: currently we have to define sim before calling DADP to avoid side effect
    nassess = 1
    sim = Simulator(pb, nassess, generation="reduction", nbins=30)
    params = District.get_sddp_solver()
    params.max_iterations = 50
    sddp = solve_SDDP(sim.model, params, 2, 1)

    pol = District.HereAndNowDP(sddp.bellmanfunctions)
    res = District.simulate(sim, pol)
elseif ALGO == "DADP"
    f, grad! = District.oracle(pb, algo)
    gdsc = @time lbfgsb(f, grad!, x0; iprint=1, pgtol=1e-5, factr=0., maxiter=40)

    pol = District.DADPPolicy([algo.models[n.name].bellmanfunctions for n in pb.nodes])
    res = District.simulate(sim, pol, )
end
