# District optimization package
# Copyright: Efficacity
# Contact: francoispacaud8@gmail.com

module District

using JLD
using JuMP
using Scenarios
using JSON

WD = pwd()


# Functions to load data from files
include("loader.jl")
# generic functions
include("model/generic.jl")
# data
include("model/data.jl")
# Definition of devices
include("model/devices.jl")

end
