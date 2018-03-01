# Macros to plot SimulationResults

NSCEN = min(100, length(res.costs))

macro ushow(res, ind, nscen=NSCEN)
    return :(plot($res.controls[:, 1:$nscen, $ind]))
end
macro xshow(res, ind, nscen=NSCEN)
    return :(plot($res.stocks[:, 1:$nscen, $ind]);)
end
macro cshow(res, nbins=50)
    return :(PyPlot.plt[:hist](.25*$res.costs', $nbins, edgecolor="k"))
end

# Plot value functions
function getgridvalues(gridb, gridh, vf)
    nx = length(gridb)
    ny = length(gridh)
    vals = zeros(Float64, nx, ny)
    for i = 1:nx, j = 1:ny
        vals[i, j] = valcuts(vf, [gridb[i], gridh[j], 20., 20.])
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
    return ax
end
