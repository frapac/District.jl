
push!(LOAD_PATH, "..")


using District, StochDynamicProgramming
using Lbfgsb, Optim

srand(2713)

ALGO = "DADP"

# Construction of the model
# Node-arc incidence matrix
A = [1. -1.]'
#= A = [-1 0; =#
#=      1 -1; =#
#=      0 1] =#


# Time span
ts = TimeSpan(180, 1)

# we build two houses
h1 = load(ts, ElecHouse(pv=4, heat=6, bat="bat0", nbins=1))
h2 = load(ts, ElecHouse(pv=0, heat=6, bat="", idhouse=2, nbins=1))
h3 = load(ts, ElecHouse(pv=4, heat=6, bat="", idhouse=2, nbins=1))
xini = Dict(h1=> [.55, 2., 20., 20.],
            h2=> [2., 20., 20.])
            #= h3=> [2., 20., 20.]) =#

# Define network
net = Network(ts, A)
# Define corresponding grid
pb = Grid(ts, [h1, h2], net)
# Build SP problems in each node
build!(pb, xini)
# Select algorithm to solve global problem
algo = DADP(pb, nsimu=1, nit=20)


# Launch Gradient Descent!
x0 = zeros(Float64, 2*95)


if ALGO == "SDDP"
    # TODO: currently we have to define sim before calling DADP to avoid side effect
    sim = Simulator(pb, 1, generation="total")
    params = District.get_sddp_solver()
    params.max_iterations = 50
    sddp = solve_SDDP(sim.model, params, 2, 1)

    pol = District.HereAndNowDP(sddp.bellmanfunctions)
    res = District.simulate(sim, pol)
elseif ALGO == "DADP"
    f, grad! = District.oracle(pb, algo)
    gdsc = @time lbfgsb(f, grad!, x0; iprint=1, pgtol=1e-5, factr=0., maxiter=40)
    pol = District.DADPPolicy([algo.models[n.name].bellmanfunctions for n in pb.nodes])
    res = District.simulate(sim, pol)
end
