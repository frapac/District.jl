
push!(LOAD_PATH, "..")

using District, StochDynamicProgramming
using Lbfgsb, Optim

#include("../../Graph.jl/generate_n_layout.jl")


srand(2713)

ALGO = "SDDP"

# Construction of the model
# Node-arc incidence matrix
#A = [1 -1]'
# graphRes = cliqueGraph(3,4)
# A = graphRes[1]
# nnodes = size(A,1)
# #= A = [-1 0; =#
# #=      1 -1; =#
# #=      0 1] =#
# #= A = [-1 0 1; =#
# #=      1 -1 0; =#
# #=      0 1 -1] =#


# Time span
ts = TimeSpan(180, 1)
# nbins = 10

# # we build two houses

# houseArray = Vector{House}(nnodes)
# houseArray[1] = load(ts, ElecHouse(pv=10, bat="", nbins=nbins, idhouse=1))
# houseArray[2] = load(ts, ElecHouse(pv=10, bat="", nbins=nbins, idhouse=2))
# houseArray[3] = load(ts, ElecHouse(pv=10, bat="", nbins=nbins, idhouse=3))

# houseArray[4] = load(ts, ElecHouse(pv=0, bat="bat0", nbins=nbins, idhouse=3))
# houseArray[5] = load(ts, ElecHouse(pv=0, bat="", nbins=nbins, idhouse=2))
# houseArray[6] = load(ts, ElecHouse(pv=4, bat="", nbins=nbins, idhouse=1))

# houseArray[7] = load(ts, ElecHouse(pv=0, bat="bat0", nbins=nbins, idhouse=1))
# houseArray[8] = load(ts, ElecHouse(pv=0, bat="", nbins=nbins, idhouse=2))
# houseArray[9] = load(ts, ElecHouse(pv=0, bat="", nbins=nbins, idhouse=3))

# houseArray[10] = load(ts, ElecHouse(pv=0, bat="", nbins=nbins, idhouse=3))
# houseArray[11] = load(ts, ElecHouse(pv=0, bat="", nbins=nbins, idhouse=1))
# houseArray[12] = load(ts, ElecHouse(pv=4, bat="bat0", nbins=nbins, idhouse=2))



# iniArray = Array{Array{Float64,1},1}(nnodes)
# for i in 1:nnodes
#     if size(houseArray[i].devices,1) >= 4
#         iniArray[i] = [.55, 2., 20., 20.]
#     else
#         iniArray[i] = [2., 20., 20.]
#     end
# end

# xini = Dict(houseArray[i]=> iniArray[i] for i in 1:nnodes)

# # Define network
# net = Network(ts, A)
# net.k2 = 1e-2 # transport cost parameters
# net.k1 = 1e-3
# # Define corresponding grid
# pb = District.Grid(ts, houseArray, net)

pb, xini = twelvehouse(nbins=10)
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
