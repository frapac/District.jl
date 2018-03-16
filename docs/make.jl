push!(LOAD_PATH, "../src")

using Documenter, District

makedocs(
    format = :html,
    sitename = "District",
    pages = [
        "Introduction" => "index.md",
        "Models" => "model.md",
        "Solvers" => "solver.md",
        "Simulation" => "simulation.md",
            ]
)
