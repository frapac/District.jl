
push!(LOAD_PATH, "..")


using District, StochDynamicProgramming
using Lbfgsb, Optim, Ipopt

srand(2713)

ALGO = "PADP"

# Construction of the model
# Node-arc incidence matrix
A = [1. -1.]'
#= A = [-1 0; =#
#=      1 -1; =#
#=      0 1] =#
#= A = [-1 0 1; =#
#=      1 -1 0; =#
#=      0 1 -1] =#


# Time span
ts = TimeSpan(180, 1)

# we build two houses
h1 = load(ts, ElecHouse(pv=4, heat=6, bat="bat0", nbins=1))
h2 = load(ts, ElecHouse(pv=0, heat=6, bat="", idhouse=2, nbins=1))
h3 = load(ts, ElecHouse(pv=4, heat=6, bat="", idhouse=2, nbins=1))
xini = Dict(h1=> [.55, 2., 20., 20.],
            h2=> [2., 20., 20.],
            h3=> [2., 20., 20.])

# Define network
net = Network(ts, A)
net.k2 = 1e-2
net.k1 = 1e-3
# Define corresponding grid
pb = Grid(ts, [h1, h2], net)
# Build SP problems in each node
build!(pb, xini, FlowInterface)
# Select algorithm to solve global problem
algo = DADP(pb, nsimu=1, nit=10)


# Launch Gradient Descent!
#= x0 = zeros(Float64, size(A, 1)*95) =#
p = EDFPrice(ts).price[1:end-1]
x0 = [p;p]


if ALGO == "SDDP"
    # TODO: currently we have to define sim before calling DADP to avoid side effect
    sim = Simulator(pb, 1, generation="total", nbins=50)
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
elseif ALGO == "IPOPT"
    f, grad! = District.oracle(pb, algo)
    nx = size(A, 1) * 95
    xL = zeros(Float64, nx)
    xU = 2*x0
    eval_g(x, g) = nothing
    eval_jac_g(x, mode, rows, cols, values ) = nothing
    prob = createProblem(nx, xL, xU, 0, Float64[], Float64[], 0, 0, f, eval_g, grad!, eval_jac_g)
    prob.x = x0
    addOption(prob, "hessian_approximation", "limited-memory")
    addOption(prob, "max_soc", 0)
    addOption(prob, "max_iter", 15)
    sol = @time solveProblem(prob)
elseif ALGO == "PADP"
    algo = PADP(pb, nsimu=1, nit=10)
    f, grad! = District.oracle(pb, algo)
    nnodes = size(A, 1)
    ntime = District.ntimes(pb) - 1
    nx = nnodes * ntime

    # x lower bound
    xL = -6*ones(Float64, nx)
    # x upper bound
    xU =  6*ones(Float64, nx)
    # constraint function
    function eval_g(x, g)
        for i in 1:ntime
            g[i] = sum(x[i:ntime:end])
        end
    end
    # gradient of constraint function
    function eval_jac_g(x, mode, rows, cols, values )
        if mode == :Structure
            for i = 1:ntime
                rows[i:ntime:end] = i
                cols[i:ntime:end] = i:ntime:(ntime*nnodes)
            end
        else
            values[:] = 1
        end
    end

    gLB = zeros(Float64, ntime)
    gUB = zeros(Float64, ntime)
    prob = createProblem(nx, xL, xU, ntime, gLB, gUB, nx, 0,
                         f, eval_g, grad!, eval_jac_g)
    prob.x = zeros(Float64, nx)
    addOption(prob, "hessian_approximation", "limited-memory")
    addOption(prob, "max_soc", 0)
    addOption(prob, "max_iter", 15)
    @time solveProblem(prob)
end
