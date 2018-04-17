
N = District.nnodes(pb)
FOLD = "$N-d"
COM = ".1000"


################################################################################
if ALGO == "SDDP"
    writecsv("results/nodal/$FOLD/sddp/exectime.csv", algo.stats.exectime)
    writecsv("results/nodal/$FOLD/sddp/lb.csv", algo.stats.lower_bounds)
    writecsv("results/nodal/$FOLD/sddp/cuts.csv", algo.bellmanfunctions)
    writecsv("results/nodal/$FOLD/sddp/incost.csv", c)

################################################################################
elseif ALGO == "DADP"

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
