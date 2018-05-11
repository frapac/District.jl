################################################################################
# District.jl
################################################################################
# Configure decomposition solvers
################################################################################
using CutPruners

"Get initial flow by solving deterministic problem"
function initialflow(pb, sddp)
    _, _, u = StochDynamicProgramming.simulate(sddp, 1)
    qind = District.getflowindex(pb)
    return u[:, 1, qind][:]
end

function runsddp(pb; nit=30)
    # TODO: currently we have to define sim before calling DADP to avoid side effect
    params = District.get_sddp_solver()
    params.max_iterations = nit
    params.compute_ub = -1
    params.reload = 20
    sddp = @time solve_SDDP(pb, params, 2, 1, prunalgo=DeMatosPruningAlgo(100))
    return sddp
end

function bfgs(pb; nsimu=1)
    algo     = DADP(pb, nsimu=nsimu, nit=30)
    f, grad! = District.oracle(pb, algo)
    p = EDFPrice(pb.ts).price[1:end-1]
    x0 = repmat(p, District.nnodes(pb))

    gdsc = @time lbfgsb(f, grad!, x0; iprint=1, pgtol=1e-5, factr=0., maxiter=50)
    return gdsc, algo
end

function qadp(pb; nsimu=1)
    algo    = QADP(pb, nsimu=nsimu, nit=20)
    f, grad! = District.oracle(pb, algo)
    narcs   = District.narcs(pb)
    ntime    = District.ntimes(pb) - 1
    nx       = narcs * ntime
    x0 = zeros(Float64, nx)

    gdsc = @time lbfgsb(f, grad!, x0; iprint=1, pgtol=1e-5, factr=0., maxiter=30)
    return gdsc, algo
end

function ipopt(pb; nsimu=1, sddpit=30, nit=10)
    algo     = DADP(pb, nsimu=nsimu, nit=sddpit)
    f, grad! = District.oracle(pb, algo)
    nx       = size(pb.net.A, 1) * 95
    xL       = zeros(Float64, nx)
    xU       = ones(Float64, nx)

    eval_g(x, g) = nothing
    eval_jac_g(x, mode, rows, cols, values ) = nothing
    prob = createProblem(nx, xL, xU, 0, Float64[], Float64[], 0, 0, f, eval_g, grad!, eval_jac_g)

    p = EDFPrice(pb.ts).price[1:end-1]
    x0 = repmat(p, District.nnodes(pb))

    addOption(prob, "hessian_approximation", "limited-memory")
    addOption(prob, "max_soc", 0)
    addOption(prob, "max_iter", nit)
    vals = Float64[]
    histmul = Float64[]
    addtrace(pb, prob, vals, histmul)

    @time solveProblem(prob)
    return prob, algo, vals, histmul
end

function quantdec(pb; nsimu=1, sddpit=20, nit=20)
    algo     = PADP(pb, nsimu=nsimu, nit=sddpit)
    f, grad! = District.oracle(pb, algo)
    nnodes   = District.nnodes(pb)
    ntime    = District.ntimes(pb) - 1
    nx       = nnodes * ntime

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

    gLB  = zeros(Float64, ntime)
    gUB  = zeros(Float64, ntime)
    prob = createProblem(nx, xL, xU, ntime, gLB, gUB, nx, 0,
                         f, eval_g, grad!, eval_jac_g)
    prob.x = zeros(Float64, nx)
    addOption(prob, "hessian_approximation", "limited-memory")
    addOption(prob, "max_soc", 0)
    addOption(prob, "max_iter", nit)
    vals = Float64[]
    histmul = Float64[]
    addtrace(pb, prob, vals, histmul)

    @time solveProblem(prob)
    return prob, algo, vals, histmul
end

function addtrace(pb, prob, vals, histmul)
    # a dirty solution to recover evolution of multipliers along iterations
    function intermediate(alg_mod::Int,
                          iter_count::Int,
                          obj_value::Float64,
                          inf_pr::Float64, inf_du::Float64,
                          mu::Float64, d_norm::Float64,
                          regularization_size::Float64,
                          alpha_du::Float64, alpha_pr::Float64,
                          ls_trials::Int)
        push!(vals, obj_value)
        push!(histmul, pb.nodes[1].conn.values...)
        return true  # Keep going
    end
    setIntermediateCallback(prob, intermediate)
end
