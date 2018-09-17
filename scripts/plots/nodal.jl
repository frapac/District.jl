using PyPlot, StatsBase, PyCall

# use seaborn to prettify plots
@pyimport seaborn
# path for results
const NNODES = "12-d"
const PATH_DATA = "results/nodal/$NNODES"
COM = ".500"


# SDDP
if false
    lb = readcsv("$PATH_DATA/sddp/lb.csv")
    tt = readcsv("$PATH_DATA/sddp/exectime.csv")
    fig = figure()
    plot(1:length(lb), lb, lw=2., c="k")
    xlabel("Iteration")
    ylabel("Lower-bound")
    seaborn.despine()
    savefig("$PATH_DATA/plots/lbsddp.pdf")

    fig = figure()
    plot(cumsum(tt), lb, lw=2., c="k")
    xlabel("Time [s]")
    ylabel("Lower-bound")
    seaborn.despine()
    savefig("$PATH_DATA/plots/lbtimesddp.pdf")
end


if false
    bw = .05
    fig = figure()
    colmap = seaborn.color_palette("Greys_r", 3)
    costs = .25*readcsv("$PATH_DATA/sddp/incosts.csv")[:, 1]
    costs_dadp = .25*readcsv("$PATH_DATA/dadp/incosts.csv")[:, 1]
    costs_qadp = .25*readcsv("$PATH_DATA/qadp/incosts.csv")[:, 1]
    cdfcost = ecdf(costs)
    xaxis = minimum(costs):.05:maximum(costs)
    #= plot(xaxis[1:end-1], diff(cdfcost.(xaxis))) =#
    seaborn.kdeplot(costs, bw=bw, lw=2., c=colmap[1], label="SDDP")
    seaborn.kdeplot(costs_dadp, bw=bw, lw=2., c=colmap[2], label="DADP")
    seaborn.kdeplot(costs_qadp, bw=bw, lw=2., c=colmap[3], label="QADP")
    seaborn.despine()
    legend()
    ylabel("Distribution")
    xlabel("Cost [€]")
    savefig("$PATH_DATA/plots/cumcost1.pdf")
end


# DADP
if false
    figure()
    λ = readcsv("$PATH_DATA/dadp$COM/1-values.csv")
    nhours = size(λ, 1)/ 4
    xtime = vec(0:.25:nhours)[1:end-1]

    haspv = [3, 5, 9, 11]
    hasbat = [1, 4, 7, 10]
    for i in 1:3
        λ = readcsv("$PATH_DATA/dadp$COM/$i-values.csv")
        if i in hasbat
            plot(xtime, λ, label="Battery", c="darkcyan")
        elseif i in haspv
            plot(xtime, λ, label="PV", c="yellowgreen")
        else
            plot(xtime, λ, label="None", c="k")
        end
    end

    for i in 4:12
        λ = readcsv("$PATH_DATA/dadp$COM/$i-values.csv")
        if i in hasbat
            plot(xtime, λ, c="darkcyan")
        elseif i in haspv
            plot(xtime, λ, c="yellowgreen")
        else
            plot(xtime, λ, c="k")
        end
    end

    plot(xtime, p[1:end-1], label="Tariff", c="coral", lw=2.5)

    xlabel("Time [h]")
    ylabel("Multipliers [€]")
    xlim(0, xtime[end])
    nhours = ceil(Int, nhours)
    xticks(0:3:nhours, 0:3:nhours)
    #= grid(which="major", axis="x") =#
    legend()
    seaborn.despine()
end


if true
    # QADP
    figure()
    λ = readcsv("$PATH_DATA/qadp$COM/1-values.csv")
    nhours = size(λ, 1)/ 4
    xtime = vec(0:.25:nhours)[1:end-1]
    haspv = [3, 5, 9, 11]
    hasbat = [1, 4, 7, 10]
    for i in 1:3
        λ = readcsv("$PATH_DATA/qadp$COM/$i-values.csv")
        if i in hasbat
            plot(xtime, λ, label="Battery", c="darkcyan")
        elseif i in haspv
            plot(xtime, λ, label="PV", c="yellowgreen")
        else
            plot(xtime, λ, label="None", c="k")
        end
    end

    for i in 4:12
        λ = readcsv("$PATH_DATA/qadp$COM/$i-values.csv")
        if i in hasbat
            plot(xtime, λ, c="darkcyan")
        elseif i in haspv
            plot(xtime, λ, c="yellowgreen")
        else
            plot(xtime, λ, c="k")
        end
    end
    xlabel("Time [h]")
    ylabel("Injection flow [kW]")
    xlim(0, xtime[end])
    nhours = ceil(Int, nhours)
    xticks(0:3:nhours, 0:3:nhours)
    #= grid(which="major", axis="x") =#
    legend()
    seaborn.despine()
    savefig("$PATH_DATA/plots/qadp-flows.pdf")
end
