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
using JSON, ProgressMeter

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
# Buildings definition
include("model/building.jl")
# Network definition
include("model/network.jl")

##########
# ALGO
include("algo/solvers.jl")
# sddp solver
include("algo/sddp.jl")

##########
# SIMULATION
# proper implementation of different policies
include("simulation/policy.jl")
# generation of assessment scenarios
include("simulation/scenarios.jl")
# Monte Carlo simulation
include("simulation/simulation.jl")

end
