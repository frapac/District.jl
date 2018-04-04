# Copyright: Efficacity
# Contact: francoispacaud8@gmail.com
################################################################################
# District.jl
################################################################################
# District optimization package
################################################################################

module District

using JLD
using JuMP
using Scenarios, StochDynamicProgramming
using JSON, ProgressMeter, MacroTools


import Base: show

WD = pwd()


# Functions to load data from files
include("loader.jl")

##########
# MODEL
# generic functions
include("model/generic.jl")
# data
include("model/data.jl")
# prices
include("model/prices.jl")
# uncertainties
include("model/uncertainties.jl")
# Definition of devices
include("model/devices.jl")
# Compute irradiance
include("model/irradiation.jl")
##########
# NODES & NETWORK
# interface between building and rest of the graph
include("model/interface.jl")
# Buildings definition
include("model/building.jl")
# Typical profile
include("model/samples.jl")
# Network definition
include("model/network.jl")
# Grid definition
include("model/grid.jl")

##########
# ALGO
include("algo/solvers.jl")
# Monte-Carlo simulation
include("algo/mcsimulation.jl")
# flow solver for graphs
include("algo/flow.jl")
# sddp solver
include("algo/sddp.jl")
# dadp solvers
include("algo/dadp.jl")
include("algo/padp.jl")
include("algo/madp.jl")
include("algo/admm.jl")

##########
# SIMULATION
# proper implementation of different policies
include("simulation/policy.jl")
# generation of assessment scenarios
include("simulation/scenarios.jl")
# Monte Carlo simulation
include("simulation/simulation.jl")


##########
# UTILS
# display functions
include("display.jl")

end
