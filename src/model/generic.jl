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

"Get starting and final timesteps of Timespan `ts`."""
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
abstract type AbstractModel end

# Definition of node in graph
abstract type AbstractNode end
