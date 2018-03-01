################################################################################
# District.jl
################################################################################
# Display functions.
################################################################################

"Print optimization model at time `t` of node `n` inside `algo`."
function dispnode(algo::AbstractDecompositionSolver, n::Int, t::Int)
    k = collect(keys(algo.models))[end:-1:1]
    println(algo.models[k[n]].solverinterface[t])
end

"Print global optimization model at time `t` inside global model."
dispglobal(sddp::SDDPInterface, t::Int) = println(sddp.solverinterface[t])

"Print global optimization model corresponding to simulation `sim`."
disppolicy(policy::AbstractPolicy, sim::Simulator, t::Int) = println(District.buildproblem!(policy, sim.model, t))
