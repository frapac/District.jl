
# use seaborn to prettify plots
# path for results
const NNODES = 3
const PATH_DATA = "results/nodal/$NNODES"


# SDDP
lb = readcsv("$PATH_DATA/sddp/lb.csv")
tt = readcsv("$PATH_DATA/sddp/exectime.csv")
@printf("SDDP time: %.2f (it: %i)\n", sum(tt)/60, length(lb))
@printf("SDDP LB: %.4f \n", lb[end] * .25)


costs = .25*readcsv("$PATH_DATA/sddp/incost.csv")[:, 1]
costs_dadp = .25*readcsv("$PATH_DATA/dadp/incosts.csv")[:, 1]
costs_qadp = .25*readcsv("$PATH_DATA/qadp/incosts.csv")[:, 1]

@printf("SDDP cost: %.3f ± %.3f\n", mean(costs), std(costs) * 1.96 / √10000)
@printf("DADP cost: %.3f ± %.3f\n", mean(costs_dadp), std(costs_dadp) * 1.96 / √10000)
@printf("QADP cost: %.3f ± %.3f\n", mean(costs_qadp), std(costs_qadp) * 1.96 / √10000)
