################################################################################
# Run Monte-Carlo estimation
################################################################################

function reload!(sddp::SDDPInterface)
    sddp.solverinterface = StochDynamicProgramming.hotstart_SDDP(sddp.spmodel,
                                         sddp.params,
                                         sddp.bellmanfunctions)
end

function tailcuts(Vs, ncuts)

    vtails = PolyhedralFunction[]

    for i in 1:endof(Vs)
        β, λ = Vs[i].betas, Vs[i].lambdas

        nc = min(ncuts, Vs[i].numCuts)

        v = PolyhedralFunction(β[1:nc], λ[1:nc, :])
        push!(vtails, v)
    end

    return vtails
end


function mc!(sddp::SDDPInterface, Vs, ncuts::Int, scen)
    # update Bellman value functions
    sddp.bellmanfunctions = tailcuts(Vs, ncuts)
    # reload JuMP Model
    reload!(sddp)

    # run Monte Carlo
    cost = StochDynamicProgramming.simulate(sddp, scen)[1]

    return mean(cost), std(cost)
end


function computeprimalMC(sddp, v, ti, to, dt, nscen=1000)
    μmc = Float64[]
    σmc = Float64[]

    srand(2713)
    scen = StochDynamicProgramming.simulate_scenarios(sddp.spmodel.noises, nscen)
    for it in ti:dt:to
        tic()
        c, s = mc!(sddp, v, it, scen)
        # we remove infinite values
        c = c[isfinite.(c)]
        tmc = toq()
        @printf("%s: %.3e", it, c)
        @printf("\t%.3e", s)
        @printf("\t%.0fs\n", tmc)
        push!(μmc, c)
        push!(σmc, s)
    end
    return μmc, σmc
end
