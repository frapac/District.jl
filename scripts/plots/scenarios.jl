################################################################################
# District.jl
################################################################################
# Scripts to plots various parameters
################################################################################

function plotdemands(scen)
    avdem = mean(scen[:, :, 2], 2)
    nhours = Int(size(scen, 1)/4)
    xtime = vec(0:.25:nhours)[1:end-1]

    fig = figure()
    plot(xtime, scen[:, :, 2], lw=.4, c="darkblue")
    plot(xtime, avdem, lw=2., c="darkred", label="Average demand")
    xlabel("Time [h]")
    ylabel("DHW demand [kW]")
    legend()
    ylim(0, 4.)
    xlim(0, xtime[end])
    xticks(0:24:nhours, 0:24:nhours)
    grid(which="major", axis="x")
end
