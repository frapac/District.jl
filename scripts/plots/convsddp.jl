using PyPlot

# path for results
const PATH_DATA = "results/nodal"

# SDDP
lb = Dict()
ct = Dict()
gap = Dict()
exectime = Dict()

for n in [2, 6, 12]
    lb[n] = readcsv("$PATH_DATA/$n/sddp/lb.csv")
    exectime[n] = readcsv("$PATH_DATA/$n/sddp/exectime.csv")
    ct2 = readcsv("$PATH_DATA/$n/sddp/algocost.csv")
    ubf2 = mean(ct2) + 1.96 * std(ct2) / √length(ct2)

    lb2 = lb[n]
    gap[n] = ubf2 ./ lb2 - 1

    @printf("Conv %i \t Gap: %.5f \t Time: %.2f\n", n, gap[n][end], sum(exectime[n])/60)
end

fig, ax = subplots()
for n in [2, 6, 12]
    plot(100*gap[n], label="$n-Nodes", lw=2.)
end
ylim(0, 2)
xlim(0, 100)
legend()
xlabel("Iteration")
ylabel("Gap [%]")
ax[:spines]["top"]["set_visible"](false)
ax[:spines]["right"]["set_visible"](false)

conv = zeros(3, 4)
for (idn, n) in enumerate([2, 6, 12])
    for (k, eps) in enumerate([0.02, 0.01, .005, .002])
        pos = findfirst(gap[n] .< eps)
        if pos > 0
            conv[idn, k] = sum(exectime[n][1:pos])
        else
            conv[idn, k] = Inf
        end
    end
end
