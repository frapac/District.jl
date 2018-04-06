using PyPlot, StatsBase, PyCall

# use seaborn to prettify plots
@pyimport seaborn
# path for results
const NNODES = 12
const PATH_DATA = "results/nodal/$NNODES"


# SDDP
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


# DADP
figure()
λ = readcsv("$PATH_DATA/dadp/1-values.csv")
nhours = size(λ, 1)/ 4
xtime = vec(0:.25:nhours)[1:end-1]
for i in 1:NNODES
    λ = readcsv("$PATH_DATA/dadp/$i-values.csv")
    plot(xtime, λ, label="Node $i")
end
xlabel("Time [h]")
ylabel("Multipliers [€]")
xlim(0, xtime[end])
xticks(0:3:nhours, 0:3:nhours)
#= grid(which="major", axis="x") =#
legend()
seaborn.despine()
savefig("$PATH_DATA/plots/dadp-multipliers.pdf")


# QADP
figure()
λ = readcsv("$PATH_DATA/qadp/1-values.csv")
nhours = size(λ, 1)/ 4
xtime = vec(0:.25:nhours)[1:end-1]
for i in 1:NNODES
    λ = readcsv("$PATH_DATA/qadp/$i-values.csv")
    plot(xtime, λ, label="Node $i")
end
xlabel("Time [h]")
ylabel("Injection flow [kW]")
xlim(0, xtime[end])
xticks(0:3:nhours, 0:3:nhours)
#= grid(which="major", axis="x") =#
legend()
seaborn.despine()
savefig("$PATH_DATA/plots/qadp-flows.pdf")
