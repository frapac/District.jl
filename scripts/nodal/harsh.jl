################################################################################
# District.jl
################################################################################
# Problems to test decomposition
# Note: Here, implement only *homogeneous* examples
################################################################################

################################################################################
function twohouse(;nbins=1)
    # Construction of the model
    tau = 1e-1
    # Node-arc incidence matrix
    A = [1. -1.]'

    # Time span
    ts = TimeSpan(180, 1)

    # we build two houses
    h1 = load(ts, ElecHouse(pv=0, heat=0, bat="bat0", nbins=nbins))
    h2 = load(ts, ElecHouse(pv=16, heat=0, bat="", idhouse=2, nbins=nbins))
    xini = Dict(h1=> [.55, 2.], h2=> [2.])

    # Define network
    net = Network(ts, A)
    net.k2 = 1e-2
    net.k1 = 1e-3
    net.maxflow[:] = 6.

    # Define corresponding grid
    return Grid(ts, [h1, h2], net), xini
end

################################################################################
function threehouse(;nbins=1)
    ts = TimeSpan(180, 1)

    h1 = load(ts, ElecHouse(pv=0, heat=0, bat="bat0", nbins=nbins))
    h2 = load(ts, ElecHouse(pv=0, heat=0, bat="", idhouse=2, nbins=nbins))
    h3 = load(ts, ElecHouse(pv=16, heat=0, bat="", idhouse=3, nbins=nbins))

    xini = Dict(h1=> [.55, 2.],
                h2=> [2.],
                h3=> [2.])

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

################################################################################
function sixhouse(;nbins=1)
    ts = TimeSpan(180, 1)

    h1 = load(ts, ElecHouse(pv=0, heat=0, bat="bat0", nbins=nbins))
    h2 = load(ts, ElecHouse(pv=0, heat=0, bat="", idhouse=2, nbins=nbins))
    h3 = load(ts, ElecHouse(pv=16, heat=0, bat="", idhouse=3, nbins=nbins))

    h4 = load(ts, ElecHouse(pv=0, heat=0, bat="bat0", nbins=nbins))
    h5 = load(ts, ElecHouse(pv=16, heat=0, bat="", idhouse=2, nbins=nbins))
    h6 = load(ts, ElecHouse(pv=0, heat=0, bat="", idhouse=3, nbins=nbins))

    xini = Dict(h1=> [.55, 2.],
                h2=> [2.],
                h3=> [2.],
                h4=> [.55, 2.],
                h5=> [2.],
                h6=> [2.])

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

################################################################################
function twelvehouse(;nbins=1)
    ts = TimeSpan(180, 1)

    h1 = load(ts, ElecHouse(pv=0, heat=0, bat="bat0", nbins=nbins))
    h2 = load(ts, ElecHouse(pv=0, heat=0, bat="", idhouse=2, nbins=nbins))
    h3 = load(ts, ElecHouse(pv=16,heat=0, bat="", idhouse=3, nbins=nbins))

    h4 = load(ts, ElecHouse(pv=0, heat=0, bat="bat0", nbins=nbins))
    h5 = load(ts, ElecHouse(pv=16,heat=0, bat="", idhouse=2, nbins=nbins))
    h6 = load(ts, ElecHouse(pv=0, heat=0, bat="", idhouse=3, nbins=nbins))

    h7 = load(ts, ElecHouse(pv=0, heat=0, bat="bat0", nbins=nbins))
    h8 = load(ts, ElecHouse(pv=0, heat=0, bat="", idhouse=2, nbins=nbins))
    h9 = load(ts, ElecHouse(pv=16,heat=0, bat="", idhouse=3, nbins=nbins))

    h10 = load(ts, ElecHouse(pv=0, heat=0, bat="bat0", nbins=nbins))
    h11 = load(ts, ElecHouse(pv=16, heat=0, bat="", idhouse=2, nbins=nbins))
    h12 = load(ts, ElecHouse(pv=0, heat=0, bat="", idhouse=3, nbins=nbins))

    xini = Dict(
                h1=> [.55, 2.],
                h2=> [2.],
                h3=> [2.],
                h4=> [.55, 2.],
                h5=> [2.],
                h6=> [2.],
                h7=> [.55, 2.],
                h8=> [2.],
                h9=> [2.],
                h10=> [.55, 2.],
                h11=> [2.],
                h12=> [2.])

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

################################################################################
function house24(;nbins=1)
    ts = TimeSpan(180, 1)

    h1 = load(ts, ElecHouse(pv=0, heat=0, bat="bat0", nbins=nbins))
    h2 = load(ts, ElecHouse(pv=0, heat=0, bat="", idhouse=2, nbins=nbins))
    h3 = load(ts, ElecHouse(pv=16, heat=0, bat="", idhouse=3, nbins=nbins))

    h4 = load(ts, ElecHouse(pv=0, heat=0, bat="bat0", nbins=nbins))
    h5 = load(ts, ElecHouse(pv=0, heat=0, bat="", idhouse=2, nbins=nbins))
    h6 = load(ts, ElecHouse(pv=16, heat=0, bat="", idhouse=3, nbins=nbins))

    h7 = load(ts, ElecHouse(pv=0, heat=0, bat="bat0", nbins=nbins))
    h8 = load(ts, ElecHouse(pv=16, heat=0, bat="", idhouse=2, nbins=nbins))
    h9 = load(ts, ElecHouse(pv=0, heat=0, bat="", idhouse=3, nbins=nbins))

    h10 = load(ts, ElecHouse(pv=0, heat=0, bat="bat0", nbins=nbins))
    h11 = load(ts, ElecHouse(pv=0, heat=0, bat="", idhouse=2, nbins=nbins))
    h12 = load(ts, ElecHouse(pv=16, heat=0, bat="", idhouse=3, nbins=nbins))

    h13 = load(ts, ElecHouse(pv=0, heat=0, bat="bat0", nbins=nbins))
    h14 = load(ts, ElecHouse(pv=16, heat=0, bat="", idhouse=2, nbins=nbins))
    h15 = load(ts, ElecHouse(pv=0, heat=0, bat="", idhouse=3, nbins=nbins))

    h16 = load(ts, ElecHouse(pv=0, heat=0, bat="bat0", nbins=nbins))
    h17 = load(ts, ElecHouse(pv=0, heat=0, bat="", idhouse=2, nbins=nbins))
    h18 = load(ts, ElecHouse(pv=16, heat=0, bat="", idhouse=3, nbins=nbins))

    h19 = load(ts, ElecHouse(pv=0, heat=0, bat="bat0", nbins=nbins))
    h20 = load(ts, ElecHouse(pv=16, heat=0, bat="", idhouse=2, nbins=nbins))
    h21 = load(ts, ElecHouse(pv=0, heat=0, bat="", idhouse=3, nbins=nbins))

    h22 = load(ts, ElecHouse(pv=0, heat=0, bat="bat0", nbins=nbins))
    h23 = load(ts, ElecHouse(pv=0, heat=0, bat="", idhouse=2, nbins=nbins))
    h24 = load(ts, ElecHouse(pv=16, heat=0, bat="", idhouse=3, nbins=nbins))
    houses = [h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11, h12,
              h13, h14, h15, h16, h17, h18, h19, h20, h21, h22, h23, h24]
    xini = Dict(
                h1=> [.55, 2.],
                h2=> [2.],
                h3=> [2.],
                h4=> [.55, 2.],
                h5=> [2.],
                h6=> [2.],
                h7=> [.55, 2.],
                h8=> [2.],
                h9=> [2.],
                h10=> [.55, 2.],
                h11=> [2.],
                h12=> [2.],
                h13=> [.55, 2.],
                h14=> [2.],
                h15=> [2.],
                h16=> [.55, 2.],
                h17=> [2.],
                h18=> [2.],
                h19=> [.55, 2.],
                h20=> [2.],
                h21=> [2.],
                h22=> [.55, 2.],
                h23=> [2.],
                h24=> [2.])


    # create matrix of other graph with LighGraph
    g = DiGraph()
    # add 24 nodes to g
    add_vertices!(g, 24)

    # add different edges for zone A
    nin = [1, 1, 1, 2, 3, 4, 4, 5, 6, 5, 7, 7, 8, 9, 10, 11]
    nout = [2, 3, 10, 3, 5, 5, 6, 6, 7, 8, 8, 9, 10, 10, 11, 12]
    for (i, j) in zip(nin, nout)
        add_edge!(g, Edge(i, j))
    end

    # add different edges for zone B
    nin = [1, 1, 2, 2, 3, 4, 5, 6, 6, 7, 9, 9, 10, 11] + 12
    nout = [2, 9, 3, 4, 5, 5, 6, 7, 8, 8, 10, 11, 11, 12] + 12
    for (i, j) in zip(nin, nout)
        add_edge!(g, Edge(i, j))
    end
    add_edge!(g, Edge(12, 13))
    add_edge!(g, Edge(2, 24))
    add_edge!(g, Edge(13, 11))
    A = -Matrix(incidence_matrix(g))

    net = Network(ts, A)
    net.k2 = 1e-2
    net.k1 = 1e-3
    net.maxflow[:] = 6.

    pb = Grid(ts, houses, net)

    return pb, xini
end


################################################################################
function house48(;nbins=1)
    ts = TimeSpan(180, 1)


    houses = []
    xini = Dict()
    for i in 1:16
        h1 = load(ts, ElecHouse(pv=0, heat=0, bat="bat0", nbins=nbins))
        h2 = load(ts, ElecHouse(pv=0, heat=0, bat="", idhouse=2, nbins=nbins))
        h3 = load(ts, ElecHouse(pv=16, heat=0, bat="", idhouse=3, nbins=nbins))

        xini[h1] = [.55, 2.]
        xini[h2] = [2.]
        xini[h3] = [2.]
        push!(houses, h1, h2, h3)
    end

    # create matrix of other graph with LighGraph
    g = DiGraph()
    # add 48 nodes to g
    add_vertices!(g, 48)

    # add different edges for zone A
    nin = [1, 1, 1, 2, 3, 4, 4, 5, 6, 5, 7, 7, 8, 9, 10, 11]
    nout = [2, 3, 10, 3, 5, 5, 6, 6, 7, 8, 8, 9, 10, 10, 11, 12]
    for (i, j) in zip(nin, nout)
        add_edge!(g, Edge(i, j))
    end

    # add different edges for zone B
    nin = [1, 1, 2, 2, 3, 4, 5, 6, 6, 7, 9, 9, 10, 11] + 12
    nout = [2, 9, 3, 4, 5, 5, 6, 7, 8, 8, 10, 11, 11, 12] + 12
    for (i, j) in zip(nin, nout)
        add_edge!(g, Edge(i, j))
    end
    add_edge!(g, Edge(12, 13))
    add_edge!(g, Edge(2, 24))
    add_edge!(g, Edge(13, 11))

    # add different edges for zone C
    nin = [1, 1, 2, 2, 3, 4, 5, 6, 6, 7, 9, 9, 10, 11] + 24
    nout = [2, 9, 3, 4, 5, 5, 6, 7, 8, 8, 10, 11, 11, 12] + 24
    for (i, j) in zip(nin, nout)
        add_edge!(g, Edge(i, j))
    end

    # add different edges for zone D
    nin = [1, 1, 1, 2, 3, 4, 4, 5, 6, 5, 7, 7, 8, 9, 10, 11] + 36
    nout = [2, 3, 10, 3, 5, 5, 6, 6, 7, 8, 8, 9, 10, 10, 11, 12] + 36
    for (i, j) in zip(nin, nout)
        add_edge!(g, Edge(i, j))
    end
    add_edge!(g, Edge(36, 37))
    add_edge!(g, Edge(26, 48))
    add_edge!(g, Edge(37, 35))

    add_edge!(g, Edge(12, 37))
    add_edge!(g, Edge(13, 25))

    A = -Matrix(incidence_matrix(g))

    net = Network(ts, A)
    net.k2 = 1e-2
    net.k1 = 0.
    net.maxflow[:] = 6.

    pb = Grid(ts, houses, net)

    return pb, xini
end
