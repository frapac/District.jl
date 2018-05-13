
N = District.nnodes(pb)
FOLD = "$N-d"
COM = ".500"


################################################################################
if ALGO == "SDDP"
    writecsv("results/nodal/$FOLD/sddp$COM/exectime.csv", algo.stats.exectime)
    writecsv("results/nodal/$FOLD/sddp$COM/lb.csv", algo.stats.lower_bounds)
    writecsv("results/nodal/$FOLD/sddp$COM/cuts.csv", algo.bellmanfunctions)
    writecsv("results/nodal/$FOLD/sddp$COM/incost.csv", c)

################################################################################
elseif ALGO == "DADP"  || ALGO == "IPOPT"

    for (inode, n) in enumerate(pb.nodes)
        sp = algo.models[n.name]
        writecsv("results/nodal/$FOLD/dadp$COM/$inode-cuts.csv", sp.bellmanfunctions)
        writecsv("results/nodal/$FOLD/dadp$COM/$inode-exectime.csv", sp.stats.exectime)
        writecsv("results/nodal/$FOLD/dadp$COM/$inode-lb.csv", sp.stats.lower_bounds)
        writecsv("results/nodal/$FOLD/dadp$COM/$inode-values.csv", n.conn.values)
    end

    writecsv("results/nodal/$FOLD/dadp$COM/algo-F.csv", algo.F)
    writecsv("results/nodal/$FOLD/dadp$COM/algo-Q.csv", algo.Q)
    writecsv("results/nodal/$FOLD/dadp$COM/net.csv", pb.net.cost)

################################################################################
elseif ALGO == "PADP" || ALGO == "QADP"

    for (inode, n) in enumerate(pb.nodes)
        sp = algo.models[n.name]
        writecsv("results/nodal/$FOLD/qadp$COM/$inode-cuts.csv", sp.bellmanfunctions)
        writecsv("results/nodal/$FOLD/qadp$COM/$inode-exectime.csv", sp.stats.exectime)
        writecsv("results/nodal/$FOLD/qadp$COM/$inode-lb.csv", sp.stats.lower_bounds)
        writecsv("results/nodal/$FOLD/qadp$COM/$inode-values.csv", n.conn.values)
    end

    writecsv("results/nodal/$FOLD/qadp$COM/algo-F.csv", algo.λ)
    writecsv("results/nodal/$FOLD/qadp$COM/algo-Q.csv", algo.μ)
    writecsv("results/nodal/$FOLD/qadp$COM/net.csv", pb.net.cost)
end
