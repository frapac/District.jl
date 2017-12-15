# Generic models

abstract type AbstractTimeSpan end

struct TimeSpan <: AbstractTimeSpan
    day::Int
    ndays::Int
end

unravel(ts::TimeSpan) = (ts.day * 96 + 1, ts.day * 96 + ts.ndays*96)
