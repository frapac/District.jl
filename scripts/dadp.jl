
push!(LOAD_PATH, "..")


using District, StochDynamicProgramming
using Lbfgsb, Optim

srand(2713)


A = [1. -1.]'


# Build problem
ts = TimeSpan(200, 1)

# we build two houses
h1 = load(ts, ElecHouse(pv=4, heat=6, bat="bat0", nbins=1))
h2 = load(ts, ElecHouse(pv=0, heat=6, bat="", idhouse=2, nbins=1))
xini = Dict(h1=> [.55, 2., 20., 20.],
            h2=> [2., 20., 20.])

net = Network(ts, A)

pb = Grid(ts, [h1, h2], net)
build!(pb, xini)
algo = DADP(pb, nsimu=1)


f, grad! = District.oracle(pb, algo)

# Launch Gradient Descent!
x0 = zeros(Float64, 190)
res = @time lbfgsb(f, grad!, x0; iprint=1, pgtol=1e-8)
#= options = Optim.Options(g_tol=1e-8, allow_f_increases=false, show_trace=true, =#
#=                         iterations=20, store_trace=true, extended_trace=false) =#
#= res = @time Optim.optimize(f, grad!, x0, BFGS(), options) =#
