
N = District.nnodes(pb)


################################################################################
if ALGO == "SDDP"
    writecsv("results/nodal/$N/sddp/exectime.csv", algo.stats.exectime)
    writecsv("results/nodal/$N/sddp/lb.csv", algo.stats.lower_bounds)
    writecsv("results/nodal/$N/sddp/cuts.csv", algo.bellmanfunctions)
    writecsv("results/nodal/$N/sddp/incost.csv", c)

################################################################################
elseif ALGO == "DADP"

    for (inode, n) in enumerate(pb.nodes)
        sp = algo.models[n.name]
        writecsv("results/nodal/$N/dadp/$inode-cuts.csv", sp.bellmanfunctions)
        writecsv("results/nodal/$N/dadp/$inode-exectime.csv", sp.stats.exectime)
        writecsv("results/nodal/$N/dadp/$inode-lb.csv", sp.stats.lower_bounds)
        writecsv("results/nodal/$N/dadp/$inode-values.csv", n.conn.values)
    end

    writecsv("results/nodal/$N/dadp/algo-F.csv", algo.F)
    writecsv("results/nodal/$N/dadp/algo-Q.csv", algo.Q)
    writecsv("results/nodal/$N/dadp/net.csv", pb.net.cost)

################################################################################
elseif ALGO == "PADP"

    for (inode, n) in enumerate(pb.nodes)
        sp = algo.models[n.name]
        writecsv("results/nodal/$N/qadp/$inode-cuts.csv", sp.bellmanfunctions)
        writecsv("results/nodal/$N/qadp/$inode-exectime.csv", sp.stats.exectime)
        writecsv("results/nodal/$N/qadp/$inode-lb.csv", sp.stats.lower_bounds)
        writecsv("results/nodal/$N/qadp/$inode-values.csv", n.conn.values)
    end

    writecsv("results/nodal/$N/qadp/algo-F.csv", algo.λ)
    writecsv("results/nodal/$N/qadp/algo-Q.csv", algo.μ)
    writecsv("results/nodal/$N/qadp/net.csv", pb.net.cost)
end
