################################################################################
# District.jl
################################################################################
# Generic models for District.jl.
# - TimeSpan is a proxy to load a specified time period.
################################################################################

export TimeSpan

################################################################################
# Define generic time period
abstract type AbstractTimeSpan end


struct TimeSpan <: AbstractTimeSpan
    # starting day
    day::Int
    # number of days considered
    ndays::Int
    # timedelta (in hours)
    δt::Float64
    # number of timesteps per day
    nt::Int
end
TimeSpan(day::Int, ndays::Int) = TimeSpan(day, ndays, .25, Int(24 / .25))

"Get starting and final timesteps of Timespan `ts`."
unravel(ts::TimeSpan) = (ts.day * ts.nt + 1, ts.day * ts.nt + ts.ndays*ts.nt)

"Get number of timesteps in TimeSpan `ts`."
ntimesteps(ts::TimeSpan) = ts.nt*ts.ndays

function weekcycle(ts::TimeSpan)
    # number of timesteps per week
    weekts = ts.nt * 7
    ti, tf = unravel(ts)
    # return modulo of ti:tf range
    return rem.((ti:tf)-1, weekts) + 1
end


################################################################################
# Define abstract type
# abstract model templates elements inside Nodes
"""
    AbstractModel

Root model for `Device` and `Uncertainties`.
A `Model` is an object intended to be an element inside a `Node`.
"""
abstract type AbstractModel end

# Definition of node in graph
"""
    AbstractNode

Generic type for `Node`.
`Node` is the basic element inside a network. `Nodes` are connected
together via edges.

# Exemple
* `House` is a node implementing a domestic house model.
"""
abstract type AbstractNode end


################################################################################
# Link between two devices
abstract type AbstractLink end
struct Link <: AbstractLink
    din::AbstractModel
    dout::AbstractModel
end

link!(n::AbstractNode, l::Link, uindex::Int, windex::Int) = link!(n, l.din, l.dout, uindex, windex)
