# Macros to plot SimulationResults
NSCEN = 10

macro ushow(res, ind, nscen=NSCEN)
    return :(plot($res.controls[:, 1:$nscen, $ind]))
end
macro xshow(res, ind, nscen=NSCEN)
    return :(plot($res.stocks[:, 1:$nscen, $ind]);)
end
macro cshow(res, nbins=50)
    return :(PyPlot.plt[:hist](.25*$res.costs', $nbins, edgecolor="k"))
end

function pplot(tab)
    fig, ax = subplots()
    plot(tab, lw=.5)
    xticks(0:12:96, 0:3:24)
    xlabel("Hour")
    ax[:spines]["top"]["set_visible"](false)
    ax[:spines]["right"]["set_visible"](false)
end

# Plot value functions
function getgridvalues(gridb, gridh, vf, to=20, ti=20)
    nx = length(gridb)
    ny = length(gridh)
    vals = zeros(Float64, nx, ny)
    for i = 1:nx, j = 1:ny
        vals[i, j] = valcuts(vf, [gridb[i], gridh[j], to, ti])
    end

    return vals
end

function valcuts(vf::PolyhedralFunction, x::Vector{Float64})
    bestcost = -Inf
    for ic in 1:vf.numCuts
        λ = @view vf.lambdas[ic, :]
        cost = dot(x, λ) + vf.betas[ic]
        if bestcost < cost
            bestcost = cost
        end
    end
    return bestcost
end

function plotvf(xx, yy, vals, ax=nothing)
    if isa(ax, Void)
        fig = figure()
        ax = fig[:add_subplot](1, 1, 1, projection="3d")
    end
    ax[:plot_surface](xx, yy, vals)
    ax[:view_init](elev=30, azim=70)
    return ax
end


function plotbellman(gridb, gridh, vf, t)
    xx = repmat(collect(gridb)', length(gridh))
    yy = repmat(collect(gridh)', length(gridb))'
    to = mean(res.stocks[t, :, 3])
    ti = mean(res.stocks[t, :, 4])

    vals = .25 * getgridvalues(gridb, gridh, vf[t], to, ti)
    plotvf(xx, yy, vals')
    xlabel("Battery level [kWh]")
    ylabel("Tank level [kWh]")
    zlabel("Value function [€]")
    tight_layout()
end

function plotcomp(xx, yy)
    fig = figure()
    ax = fig[:add_subplot](1, 1, 1, projection="3d")
    ax[:plot_surface](xx, yy, .25vals_p', label="Price function")
    ax[:plot_surface](xx, yy, .25vals_q', label="Quant function")
    xlabel("Battery level [kWh]")
    ylabel("Tank level [kWh]")
    zlabel("Value function [€]")
end
