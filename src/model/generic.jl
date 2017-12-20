# Generic models
# TODO: change files' name

export TimeSpan

abstract type AbstractTimeSpan end

struct TimeSpan <: AbstractTimeSpan
    day::Int
    ndays::Int
    δt::Float64
    nt::Int
end
TimeSpan(day, ndays) = TimeSpan(day, ndays, .25, Int(24 / .25))

unravel(ts::TimeSpan) = (ts.day * ts.nt + 1, ts.day * ts.nt + ts.ndays*ts.nt)
ntimesteps(ts::TimeSpan) = ts.nt*ts.ndays

function weekcycle(ts::TimeSpan)
    # number of timesteps per week
    weekts = ts.nt * 7
    ti, tf = unravel(ts)
    return rem.(ti:tf, weekts)
end
