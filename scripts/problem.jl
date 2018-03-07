
function twohouse(;nbins=1)
    # Construction of the model
    tau = 1e-1
    # Node-arc incidence matrix
    A = [1. -1.]'


    # Time span
    ts = TimeSpan(180, 1)

    # we build two houses
    h1 = load(ts, ElecHouse(pv=4, heat=6, bat="bat0", nbins=nbins))
    h2 = load(ts, ElecHouse(pv=0, heat=6, bat="", idhouse=2, nbins=nbins))
    xini = Dict(h1=> [.55, 2., 20., 20.], h2=> [2., 20., 20.])

    # Define network
    net = Network(ts, A)
    net.k2 = 1e-2
    net.k1 = 1e-3
    net.maxflow[:] = 6.

    # Define corresponding grid
    return Grid(ts, [h1, h2], net), xini
end


function threehouse(;nbins=1)
    ts = TimeSpan(180, 1)

    h1 = load(ts, ElecHouse(pv=4, heat=6, bat="bat0", nbins=nbins))
    h2 = load(ts, ElecHouse(pv=0, heat=6, bat="", idhouse=2, nbins=nbins))
    h3 = load(ts, ElecHouse(pv=4, heat=6, bat="", idhouse=3, nbins=nbins))

    xini = Dict(h1=> [.55, 2., 20., 20.],
                h2=> [2., 20., 20.],
                h3=> [2., 20., 20.])

    A = [-1 0 1;
         1 -1 0;
         0 1 -1]
    net = Network(ts, A)
    net.k2 = 1e-2
    net.k1 = 1e-3
    net.maxflow[:] = 6.

    pb = Grid(ts, [h1, h2, h3], net)

    return pb, xini
end

function sixhouse(;nbins=1)
    ts = TimeSpan(180, 1)

    h1 = load(ts, ElecHouse(pv=4, heat=6, bat="bat0", nbins=nbins))
    h2 = load(ts, ElecHouse(pv=0, heat=6, bat="", idhouse=2, nbins=nbins))
    h3 = load(ts, ElecHouse(pv=4, heat=6, bat="", idhouse=3, nbins=nbins))

    h4 = load(ts, ElecHouse(pv=4, heat=6, bat="bat0", nbins=nbins))
    h5 = load(ts, ElecHouse(pv=0, heat=6, bat="", idhouse=2, nbins=nbins))
    h6 = load(ts, ElecHouse(pv=4, heat=6, bat="", idhouse=3, nbins=nbins))

    xini = Dict(h1=> [.55, 2., 20., 20.],
                h2=> [2., 20., 20.],
                h3=> [2., 20., 20.],
                h4=> [.55, 2., 20., 20.],
                h5=> [2., 20., 20.],
                h6=> [2., 20., 20.])

    A = [1 0 -1 0 0 0 0;
         -1 1 0 0 0 0 0;
         0 -1 1 0 0 0 1;
         0  0 0 1 -1 0 0;
         0  0 0 -1 0 1 -1;
         0  0 0 0  1 -1 0]
    net = Network(ts, A)
    net.k2 = 1e-2
    net.k1 = 1e-3
    net.maxflow[:] = 6.

    pb = Grid(ts, [h1, h2, h3, h4, h5, h6], net)

    return pb, xini
end

function twelvehouse(;nbins=1)
    ts = TimeSpan(180, 1)

    h1 = load(ts, ElecHouse(pv=4, heat=6, bat="bat0", nbins=nbins))
    h2 = load(ts, ElecHouse(pv=0, heat=6, bat="", idhouse=2, nbins=nbins))
    h3 = load(ts, ElecHouse(pv=4, heat=6, bat="", idhouse=3, nbins=nbins))

    h4 = load(ts, ElecHouse(pv=4, heat=6, bat="bat0", nbins=nbins))
    h5 = load(ts, ElecHouse(pv=0, heat=6, bat="", idhouse=2, nbins=nbins))
    h6 = load(ts, ElecHouse(pv=4, heat=6, bat="", idhouse=3, nbins=nbins))

    h7 = load(ts, ElecHouse(pv=4, heat=6, bat="bat0", nbins=nbins))
    h8 = load(ts, ElecHouse(pv=0, heat=6, bat="", idhouse=2, nbins=nbins))
    h9 = load(ts, ElecHouse(pv=4, heat=6, bat="", idhouse=3, nbins=nbins))

    h10 = load(ts, ElecHouse(pv=4, heat=6, bat="bat0", nbins=nbins))
    h11 = load(ts, ElecHouse(pv=0, heat=6, bat="", idhouse=2, nbins=nbins))
    h12 = load(ts, ElecHouse(pv=4, heat=6, bat="", idhouse=3, nbins=nbins))

    xini = Dict(
                h1=> [.55, 2., 20., 20.],
                h2=> [2., 20., 20.],
                h3=> [2., 20., 20.],
                h4=> [.55, 2., 20., 20.],
                h5=> [2., 20., 20.],
                h6=> [2., 20., 20.],
                h7=> [.55, 2., 20., 20.],
                h8=> [2., 20., 20.],
                h9=> [2., 20., 20.],
                h10=> [.55, 2., 20., 20.],
                h11=> [2., 20., 20.],
                h12=> [2., 20., 20.])

    A1 = Float64[1 0 -1 0 0 0 0;
         -1 1 0 0 0 0 0;
         0 -1 1 0 0 0 1;
         0  0 0 1 -1 0 0;
         0  0 0 -1 0 1 -1;
         0  0 0 0  1 -1 0]
    A2 = Float64[1 1 0 0 0 0;
          -1 0 1 0 0 0;
          0 -1 0 1 0 0;
          0 0 -1 -1 1 0;
          0 0 0 0 -1 1;
          0 0 0 0 0 -1]
    O1 = zeros(Float64, 6, 6)
    O2 = zeros(Float64, 6, 7)

    A = [A1 O1; O2 A2]
    A = [A zeros(Float64, 12, 3)]
    A[1, 14] = 1
    A[10, 14] = -1
    A[5, 15] = 1
    A[8, 15] = -1
    A[6, 16] = 1
    A[7, 16] = -1

    net = Network(ts, A)
    net.k2 = 1e-2
    net.k1 = 1e-3
    net.maxflow[:] = 6.

    pb = Grid(ts, [h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11, h12], net)

    return pb, xini
end
