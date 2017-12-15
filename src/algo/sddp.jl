

""" Build LP model with pre-computed cuts.

$(SIGNATURES)

# Arguments
* `model::SPModel`
* `param::SDDPparameters`
* `V::Vector{PolyhedralFunction}`
* `relaxation::Bool`
* `form::string`
  Specify how to approximate the online expectancy
"""
function build_models(model, param::SDDPparameters,
                      V::Vector{StochDynamicProgramming.PolyhedralFunction};
                      relaxation=true, modelling=:pred)

    models = Vector{JuMP.Model}(model.stageNumber-1)
    p = Progress(model.stageNumber-1, 1)

    for t = 1:model.stageNumber-1
        if modelling == :pred
            models[t] = build_model(model, param, t, V, relaxation)
        else
            models[t] = build_model_dh(model, param, t, V)
        end
        next!(p)
    end
    try
        model.finalCost(model, models[end])
    catch
        println("WARN: fail to define final cost")
    end
    return models
end


# TODO: fix final costs in definition of costs
function build_model(model::StochDynamicProgramming.SPModel,
                     param::StochDynamicProgramming.SDDPparameters,
                     t::Int64,
                     V::Vector{StochDynamicProgramming.PolyhedralFunction}, relaxation::Bool)
    m = Model(solver=param.SOLVER)

    nx = model.dimStates
    nu = model.dimControls
    nw = model.dimNoises

    @variable(m,  model.xlim[i][1] <= x[i=1:nx] <= model.xlim[i][2])
    @variable(m,  model.xlim[i][1] <= xf[i=1:nx]<= model.xlim[i][2])
    @variable(m,  model.ulim[i][1] <= u[i=1:nu] <=  model.ulim[i][2],
              category=model.controlCat[i])
    @variable(m, alpha)

    @variable(m, w[1:nw] == 0)
    m.ext[:cons] = @constraint(m, state_constraint, x .== 0)

    @constraint(m, xf .== model.dynamics(t, x, u, w))

    if isa(model.costFunctions, Function)
    try
        @objective(m, Min, model.costFunctions(t, x, u, w) + alpha)
    catch
        @objective(m, Min, model.costFunctions(m, t, x, u, w) + alpha)
    end

    elseif isa(model.costFunctions, Vector{Function})
        @variable(m, cost)

        for i in 1:length(model.costFunctions)
            @constraint(m, cost >= model.costFunctions[i](t, x, u, w))
        end
        @objective(m, Min, cost + alpha)
    end

    for nc in 1:V[t+1].numCuts
        lambda = vec(V[t+1].lambdas[nc, :])
        @constraint(m, V[t+1].betas[nc] + dot(lambda, xf) <= alpha)
    end
    return m
end



function build_model_dh(model, param, t, V)
    m = Model(solver=param.SOLVER)
    law = model.noises

    nx = model.dimStates
    nu = model.dimControls
    nw = model.dimNoises

    ns = law[t].supportSize
    ξ = collect(law[t].support[:, :])
    πp = law[t].proba

    @variable(m, model.xlim[i][1] <= x[i=1:nx] <= model.xlim[i][2])
    @variable(m, model.ulim[i][1] <= u[i=1:nu] <=  model.ulim[i][2])
    @variable(m, model.xlim[i][1] <= xf[i=1:nx, j=1:ns]<= model.xlim[i][2])
    @variable(m, alpha[1:ns])

    @variable(m, w[1:nw] == 0)
    m.ext[:cons] = @constraint(m, state_constraint, x .== 0)

    for j=1:ns
        @constraint(m, xf[:, j] .== model.dynamics(t, x, u, ξ[:, j]))
    end

    # add objective as minimization of expectancy:
    @objective(m, Min,
               sum(πp[j]*(model.costFunctions(m, t, x, u, ξ[:, j]) +
                    alpha[j]) for j in 1:ns))

    for nc in 1:V[t+1].numCuts
        lambda = vec(V[t+1].lambdas[nc, :])
        for j=1:ns
            @constraint(m, V[t+1].betas[nc] + dot(lambda, xf[:, j]) <= alpha[j])
        end
    end

    return m
end


