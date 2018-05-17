using PyPlot

# path for results
COM = ".100"
const PATH_DATA = "results/nodal"

# SDDP
lb = Dict()
ct = Dict()
gap = Dict()
exectime = Dict()

for n in [3, 6, 12, 24]
    lb[n] = readcsv("$PATH_DATA/$n-d/sddp$COM/lb.csv")
    exectime[n] = readcsv("$PATH_DATA/$n-d/sddp$COM/exectime.csv")
    ct2 = readcsv("$PATH_DATA/$n-d/sddp$COM/incost.csv")
    ubf2 = mean(ct2) + 1.96 * std(ct2) / √length(ct2)

    lb2 = lb[n]
    gap[n] = ubf2 ./ lb2 - 1

    @printf("Conv %i \t Gap: %.5f \t Time: %.2f\n", n, gap[n][end], sum(exectime[n])/60)
end

fig, ax = subplots()
for n in [3, 6, 12, 24]
    plot(100*gap[n], label="$n-Nodes", lw=2.)
end
ylim(0, 2)
xlim(0, 100)
legend()
xlabel("Iteration")
ylabel("Gap [%]")
ax[:spines]["top"]["set_visible"](false)
ax[:spines]["right"]["set_visible"](false)

conv = zeros(4, 4)
iter = zeros(4, 4)
for (idn, n) in enumerate([3, 6, 12, 24])
    for (k, eps) in enumerate([0.02, 0.01, .005, .002])
        pos = findfirst(gap[n] .< eps)
        if pos > 0
            conv[idn, k] = sum(exectime[n][1:pos])
            iter[idn, k] = pos
        else
            conv[idn, k] = Inf
            iter[idn, k] = Inf
        end
    end
end
