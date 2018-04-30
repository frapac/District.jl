
"""Display convergence of SDDP in primal and in dual."""
function dispconv(ubp, stdp; nscen=1000, ptol=.95, delta=10)
    ΔI= delta
    MAXIT = sddp.stats.niterations
    tol = √2 * erfinv(2*ptol - 1)
    tol2 = √2 * erfinv(2*.999 - 1)
    ubf = sddp.stats.upper_bounds
    lprim = .25sddp.stats.lower_bounds

    fig, ax = subplots(figsize=(9,6))

    plot(lprim, c="darkred", lw=4, label="SDDP LB")
    #= plot(ubf, c="k", lw=.1, label="Forward primal cost") =#


    plot(ΔI:ΔI:MAXIT, ubp, color="k", lw=1.5, label="SDDP UB", marker="s")
    plot(ΔI:ΔI:MAXIT, ubp + tol*stdp/√nscen, color="k", lw=1, linestyle="--")
    plot(ΔI:ΔI:MAXIT, ubp - tol*stdp/√nscen, color="k", lw=1, linestyle="--")
    #= plot(ΔI:ΔI:MAXIT, ubp + tol2*stdp/√nscen, color="k", lw=1, linestyle=":") =#
    #= plot(ΔI:ΔI:MAXIT, ubp - tol2*stdp/√nscen, color="k", lw=1, linestyle=":") =#
    fill_between(ΔI:ΔI:MAXIT, ubp - tol*stdp/√nscen, ubp + tol*stdp/√nscen,
                 color="grey", alpha=.4, label="Confidence ($(100*ptol)%)")
    #= fill_between(ΔI:ΔI:MAXIT, ubp - tol2*stdp/√nscen, ubp + tol2*stdp/√nscen, =#
    #=              color="grey", alpha=.1, label="Confidence ($(100*.999)%)") =#


    ax[:spines]["top"]["set_visible"](false)
    ax[:spines]["right"]["set_visible"](false)
    legend(loc=4)
    xlabel("Iterations")
    ylabel("Cost [€]")
    tight_layout()
end


function disptrajectories(res1, res2, nx; nscen=100)
    fig = figure(figsize=(10, 5))
    ax = subplot(1, 2, 1)
    plot(res1.stocks[:, 1:nscen, nx], c="darkslateblue", lw=.7);
    xlabel("Hours", fontsize=15)
    ylabel("Charge [kWh]", fontsize=15)
    tick_params(labelsize=14)
    title("SDDP")

    ax2 = subplot(1, 2, 2, sharex=ax,  sharey=ax)
    plot(res2.stocks[:, 1:nscen, nx], c="darkslateblue", lw=.7);

    xticks(0:12:96, 0:3:24)
    xlim(0, 96)
    xlabel("Hours", fontsize=15)
    tight_layout()
    tick_params(labelsize=14)
    title("MPC")
end


function dispub(lb, ubp, stdp; nscen=1000, ptol=.95, delta=10, nit=300)
    ΔI= delta
    MAXIT = nit
    tol = √2 * erfinv(2*ptol - 1)
    tol2 = √2 * erfinv(2*.999 - 1)

    fig, ax = subplots(figsize=(9,6))
    plot(lb, c="darkred", lw=4, label="SDDP LB")


    plot(ΔI:ΔI:MAXIT, ubp, color="k", lw=1.5, label="SDDP UB", marker="s")
    plot(ΔI:ΔI:MAXIT, ubp + tol*stdp/√nscen, color="k", lw=1, linestyle="--")
    plot(ΔI:ΔI:MAXIT, ubp - tol*stdp/√nscen, color="k", lw=1, linestyle="--")
    #= plot(ΔI:ΔI:MAXIT, ubp + tol2*stdp/√nscen, color="k", lw=1, linestyle=":") =#
    #= plot(ΔI:ΔI:MAXIT, ubp - tol2*stdp/√nscen, color="k", lw=1, linestyle=":") =#
    fill_between(ΔI:ΔI:MAXIT, ubp - tol*stdp/√nscen, ubp + tol*stdp/√nscen,
                 color="grey", alpha=.4, label="Confidence ($(100*ptol)%)")
    #= fill_between(ΔI:ΔI:MAXIT, ubp - tol2*stdp/√nscen, ubp + tol2*stdp/√nscen, =#
    #=              color="grey", alpha=.1, label="Confidence ($(100*.999)%)") =#


    ax[:spines]["top"]["set_visible"](false)
    ax[:spines]["right"]["set_visible"](false)
    legend(loc=4)
    xlabel("Iterations")
    ylabel("Cost [€]")
    tight_layout()
end
