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

function plotweather(ts)
    nhours = ts.ndays * 24
    xtime = vec(0:.25:nhours)[1:end-1]

    te = loadweather(OutdoorTemperature(), ts)
    ghi = loadweather(District.GHI(), ts)
    ghi_cs = 4 * District.loaddata(ts, 14)

    fig, ax = subplots(nrows=2, ncols=1, sharex=true, figsize=(10,6))

    ax[1][:plot](xtime, te, lw=.7, c="k")
    ax[1][:set_ylabel]("Outdoor temperature [°C]")

    ax[2][:plot](xtime, ghi / 1000, lw=.7, c="k", label="GHI")
    ax[2][:plot](xtime, ghi_cs / 1000, lw=.7, c="r", linestyle="--", label="Clear Sky GHI")
    ax[2][:set_ylabel]("Radiation [kW / m2]")
    ax[2][:legend](loc="upper right")

    xlim(0, xtime[end])
    xticks(0:24:nhours, 0:24:nhours)
    ax[2][:grid](which="major", axis="x", linestyle=":", lw=.5)
    ax[1][:grid](which="major", axis="x", linestyle=":", lw=.5)
    xlabel("Time [h]")
    tight_layout()
end
