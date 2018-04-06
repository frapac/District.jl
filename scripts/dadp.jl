
push!(LOAD_PATH, "..")


using District, StochDynamicProgramming
using Lbfgsb, Optim, Ipopt

srand(2713)

ALGO = "DADP"

include("graph.jl")
include("problem.jl")

# Construction of the model
# Node-arc incidence matrix
# Time span
ts = TimeSpan(180, 1)

# Noise discretization
nbins = 1

# Build grid
pb, xini = twelvehouse(nbins=nbins)

nnodes = size(pb.nodes,1)

# Build SP problems in each node
build!(pb, xini, PriceInterface)


# Select algorithm to solve global problem
algo = DADP(pb, nsimu=1, nit=10)


# Launch Gradient Descent!
#= x0 = zeros(Float64, size(A, 1)*95) =#
p = EDFPrice(ts).price[1:end-1]
mul0 = p
for i in 1:nnodes-1
    mul0 = vcat(mul0, p)
end


if ALGO == "SDDP"
    nassess = 1
    # TODO: currently we have to define sim before calling DADP to avoid side effect
    sim = Simulator(pb, nassess, generation="reduction", nbins=1)
    params = District.get_sddp_solver()
    params.max_iterations = 5
    sddp = solve_SDDP(sim.model, params, 2, 1)

    pol = District.HereAndNowDP(sddp.bellmanfunctions)
    res = District.simulate(sim, pol)
elseif ALGO == "DADP"
    f, grad! = District.oracle(pb, algo)
    gdsc = @time lbfgsb(f, grad!, mul0; iprint=1, pgtol=1e-5, factr=0., maxiter=10)
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
