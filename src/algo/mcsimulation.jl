################################################################################
# District.jl
################################################################################
# Extract sensitivity w.r.t. a given constraint.
# Estimation by Monte Carlo sampling.
################################################################################

############################################################
# Primal estimation
"""
    mcsimulation(sddp::StochDynamicProgramming.SDDPInterface, scenarios::Array{Float64, 3})

Compute average importation in node by Monte-Carlo, along scenarios
`scenarios`.

# Returns
* `costs::Vector{Float64}` size=(nscen,)
Costs along assessment scenarios.

* `import::Array{Float64, 2}` size=(ntime, nscen)
Importation along assessment scenarios.
"""
function mcsimulation(sddp::StochDynamicProgramming.SDDPInterface, scenarios::Array{Float64, 3}, ninjection::Int64)
    return mcsimulation(sddp.spmodel, sddp.params, sddp.solverinterface, scenarios, ninjection)
end

function mcsimulation(model::StochDynamicProgramming.SPModel,
                      param::StochDynamicProgramming.SDDPparameters,
                      solverProblems::Vector{JuMP.Model},
                      scenarios::Array{Float64},
                      ninj::Int64)
    

    T = model.stageNumber
    nb_forward = size(scenarios, 2)

    stockTrajectories = zeros(T, nb_forward, model.dimStates)
    # import is last position in controls
    importation = zeros(ninj, T - 1, nb_forward)

    # Set first value of stocks equal to x0:
    for k in 1:nb_forward
        stockTrajectories[1, k, :] = model.initialState
    end

    # Store costs of different scenarios in an array:
    costs = zeros(nb_forward)

    for t=1:T-1
        for k = 1:nb_forward
            # Collect current state and noise:
            xt = stockTrajectories[t, k, :]
            ξt = scenarios[t, k, :]
            
            sol, ts = StochDynamicProgramming.solve_one_step_one_alea(model, param,
                                              solverProblems[t], t, xt, ξt)

            if sol.status
                # extract sensitivity
                importation[:, t, k] = sol.uopt[end-ninj+1:end]

                stockTrajectories[t+1, k, :] = sol.xf
                # and the current cost:
                costs[k] += sol.objval - sol.θ
                if t==T-1
                    costs[k] += sol.θ
                end       
            else
                costs[k] = Inf
                break
            end
        end
    end
    return costs, importation
end



############################################################
# Dual estimation
"""
    qsensitivity(sddp::SDDPInterface, scenarios::Array{Float64, 3}

Extract sensitivity of problem w.r.t. coupling constraints along `scenarios`.
Estimation by Monte-Carlo sampling.
"""
function qsensitivity(sddp::StochDynamicProgramming.SDDPInterface, scenarios::Array{Float64, 3})
    return sensitivity(sddp.spmodel, sddp.params, sddp.solverinterface,
                       scenarios, :coupling)
end

function sensitivity(model::StochDynamicProgramming.SPModel,
                     param::StochDynamicProgramming.SDDPparameters,
                     solverProblems::Vector{JuMP.Model},
                     scenarios::Array{Float64},
                     refcons::Symbol=:cons)


    T = model.stageNumber
    nb_forward = size(scenarios, 2)


    stockTrajectories = zeros(T, nb_forward, model.dimStates)

    cons = solverProblems[1].ext[refcons]
    np = isa(cons, SimpleVector) ? length(cons) : 1
    sensitivity = zeros(T - 1, nb_forward, np)

    # Set first value of stocks equal to x0:
    for k in 1:nb_forward
        stockTrajectories[1, k, :] = model.initialState
    end

    # Store costs of different scenarios in an array:
    costs = zeros(nb_forward)

    for t=1:T-1
        for k = 1:nb_forward
            # Collect current state and noise:
            xt = stockTrajectories[t, k, :]
            ξt = scenarios[t, k, :]

            sol, ts = StochDynamicProgramming.solve_one_step_one_alea(model, param,
                                              solverProblems[t], t, xt, ξt)


            if sol.status
                # extract sensitivity
                sensitivity[t, k, :] = JuMP.getdual(solverProblems[t].ext[refcons])
                stockTrajectories[t+1, k, :] = sol.xf
                # and the current cost:
                costs[k] += sol.objval - sol.θ
                if t==T-1
                    costs[k] += sol.θ
                end
            else
                costs[k] = Inf
                break
            end
        end
    end

    # remove unvalid trajectories
    isvalid = isfinite.(costs)
    return costs[isvalid], sensitivity[:, isvalid]
end
