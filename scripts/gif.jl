
for i in 1:30
    fig, ax = subplots()
    plot(l[:, i], lw=1., c="k")
    xticks(0:12:96, 0:3:24)
    xlabel("Hour")
    ylim(0.12, 0.165)
    ax[:spines]["top"]["set_visible"](false)
    ax[:spines]["right"]["set_visible"](false)
    savefig("results/pictures/gif/mul$i.jpg")
    close()
end
