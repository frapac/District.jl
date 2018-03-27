################################################################################
# District.jl
################################################################################
# Scripts to plots simulation results
################################################################################

push!(LOAD_PATH, "..")

using PyPlot, JLD, District

NSCEN = 60
include("../problem.jl")
include("graph.jl")


pb, xini = twohouse(nbins=10)

N = District.nnodes(pb)

πel = EDFPrice(pb.ts).price

# get correspondence to label plots
xname, uname, wname = District.getcorrespondance(pb)

# load database storing results
db = load("results/nodal/$N/simuls.jld")
res = db["res"]
eld = db["load"]


ntime, nsimu, nx = size(res.stocks)
nhours = Int(ntime/4)
xtime = vec(0:.25:nhours)[1:end-1]

################################################################################
# Plot stocks level
nx = 1
for i in 2:nx
    series = res.stocks[:, :, i]
    avdem = mean(series, 2)

    fig, ax = subplots()
    plot(xtime, series[:, 1:NSCEN], lw=.4, c="darkblue")
    plot(xtime, avdem, lw=6., c="k")
    xlabel("Time [h]")
    ylabel(xname[i])
    xlim(0, xtime[end])
    xticks(0:3:nhours, 0:3:nhours)

    if contains(xname[i], "Battery")
        ylim(.5, 3.)
    end
    ax[:spines]["top"]["set_visible"](false)
    ax[:spines]["right"]["set_visible"](false)
end

################################################################################
# Plot importation
series = .25 * eld
avdem = mean(series, 2)
fig, ax = subplots()
plot(xtime, series[:, 1:NSCEN], lw=.4, c="darkblue")
plot(xtime, avdem, lw=6., c="k")
xlabel("Time [h]")
ylabel("Elec import. [kWh]")
xlim(0, xtime[end])
xticks(0:3:nhours, 0:3:nhours)

ax[:spines]["top"]["set_visible"](false)
ax[:spines]["right"]["set_visible"](false)
savefig("/home/fpacaud/Documents/RENDU/RAPPORTS/DADP/img/assess-$N-import.pdf")


################################################################################
# Display graph and exchanges
na = District.narcs(pb)
q = mean(res.controls[:, :, end-na+1:end], 2)
plotflow(pb.net.A, mean(q, 1), darrow=true)
savefig("/home/fpacaud/Documents/RENDU/RAPPORTS/DADP/img/assess-$N-graph.pdf")

################################################################################
elcost = sum(mean(πel .* eld, 2)[1:end-1])
@printf("Costs\n")
@printf("Total cost: \t %.2f €\n", .25*mean(res.costs))
@printf("Elec cost: \t %.2f €\n", .25*elcost)
@printf("Import: \t %.2f kWh\n", .25*sum(mean(eld, 2)))
