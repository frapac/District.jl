# Copyright Ian Dunning
# https://github.com/IainNZ/GraphLayout.jl/blob/master/src/stress.jl

"""
Compute graph layout using stress majorization

Inputs:

    δ: Matrix of pairwise distances
    p: Dimension of embedding (default: 2)
    w: Matrix of weights. If not specified, defaults to
           w[i,j] = δ[i,j]^-2 if δ[i,j] is nonzero, or 0 otherwise
    X0: Initial guess for the layout. Coordinates are given in rows.
        If not specified, default to random matrix of Gaussians

Additional optional keyword arguments control the convergence of the algorithm
and the additional output as requested:

    maxiter:   Maximum number of iterations. Default: 400size(X0, 1)^2
    abstols:   Absolute tolerance for convergence of stress.
               The iterations terminate if the difference between two
               successive stresses is less than abstol.
               Default: √(eps(eltype(X0))
    reltols:   Relative tolerance for convergence of stress.
               The iterations terminate if the difference between two
               successive stresses relative to the current stress is less than
               reltol. Default: √(eps(eltype(X0))
    abstolx:   Absolute tolerance for convergence of layout.
               The iterations terminate if the Frobenius norm of two successive
               layouts is less than abstolx. Default: √(eps(eltype(X0))
    verbose:   If true, prints convergence information at each iteration.
               Default: false
    returnall: If true, returns all iterates and their associated stresses.
               If false (default), returns the last iterate

Output:

    The final layout X, with coordinates given in rows, unless returnall=true.

Reference:

    The main equation to solve is (8) of:

    @incollection{
        author = {Emden R Gansner and Yehuda Koren and Stephen North},
        title = {Graph Drawing by Stress Majorization}
        year={2005},
        isbn={978-3-540-24528-5},
        booktitle={Graph Drawing},
        seriesvolume={3383},
        series={Lecture Notes in Computer Science},
        editor={Pach, J\'anos},
        doi={10.1007/978-3-540-31843-9_25},
        publisher={Springer Berlin Heidelberg},
        pages={239--250},
    }
"""
function layout_stressmajorize_adj(δ, p::Int=2, w=nothing, X0=randn(size(δ, 1), p);
        maxiter = 400size(X0, 1)^2, abstols=√(eps(eltype(X0))),
        reltols=√(eps(eltype(X0))), abstolx=√(eps(eltype(X0))),
        verbose = false, returnall = false)

    @assert size(X0, 2)==p

    if w==nothing
        w = δ.^-2
        w[!isfinite(w)] = 0
    end

    @assert size(X0, 1)==size(δ, 1)==size(δ, 2)==size(w, 1)==size(w, 2)
    Lw = weightedlaplacian(w)
    pinvLw = pinv(Lw)
    newstress = stress(X0, δ, w)
    Xs = Matrix[X0]
    stresses = [newstress]
    iter = 0
    for iter = 1:maxiter
        #TODO the faster way is to drop the first row and col from the iteration
        X = pinvLw * (LZ(X0, δ, w)*X0)
        @assert all(isfinite(X))
        newstress, oldstress = stress(X, δ, w), newstress
        verbose && info("""Iteration $iter
        Change in coordinates: $(vecnorm(X - X0))
        Stress: $newstress (change: $(newstress-oldstress))
        """)
        push!(Xs, X)
        push!(stresses, newstress)
        abs(newstress - oldstress) < reltols * newstress && break
        abs(newstress - oldstress) < abstols && break
        vecnorm(X - X0) < abstolx && break
        X0 = X
    end
    iter == maxiter && warn("Maximum number of iterations reached without convergence")
    returnall ? (Xs, stresses) : Xs[end]
end

"""
Stress function to majorize

Input:
    X: A particular layout (coordinates in rows)
    d: Matrix of pairwise distances
    w: Weights for each pairwise distance

See (1) of Reference
"""
function stress(X, d=fill(1.0, size(X, 1), size(X, 1)), w=nothing)
    s = 0.0
    n = size(X, 1)
    if w==nothing
        w = d.^-2
        w[!isfinite(w)] = 0
    end
    @assert n==size(d, 1)==size(d, 2)==size(w, 1)==size(w, 2)
    for j=1:n, i=1:j-1
        s += w[i, j] * (norm(X[i,:] - X[j,:]) - d[i,j])^2
    end
    @assert isfinite(s)
    s
end

"""
Compute weighted Laplacian given ideal weights w

Lʷ defined in (4) of the Reference
"""
function weightedlaplacian(w)
    n = Base.LinAlg.chksquare(w)
    T = eltype(w)
    Lw = zeros(T, n, n)
    for i=1:n
        D = zero(T)
        for j=1:n
            i==j && continue
            Lw[i, j] = -w[i, j]
            D += w[i, j]
        end
        Lw[i, i]=D
    end
    Lw
end

"""
Computes L^Z defined in (5) of the Reference

Input: Z: current layout (coordinates)
       d: Ideal distances (default: all 1)
       w: weights (default: d.^-2)
"""
function LZ(Z, d, w)
    n = size(Z, 1)
    L = zeros(n, n)
    for i=1:n
        D = 0.0
        for j=1:n
            i==j && continue
            nrmz = norm(Z[i,:] - Z[j,:])
            nrmz==0 && continue
            δ = w[i, j] * d[i, j]
            L[i, j] = -δ/nrmz
            D -= -δ/nrmz
        end
        L[i, i] = D
    end
    @assert all(isfinite(L))
    L
end

"""
    Use the spring/repulsion model of Fruchterman and Reingold (1991):
        Attractive force:  f_a(d) =  d^2 / k
        Repulsive force:  f_r(d) = -k^2 / d
    where d is distance between two vertices and the optimal distance
    between vertices k is defined as C * sqrt( area / num_vertices )
    where C is a parameter we can adjust

    Arguments:
    adj_matrix Adjacency matrix of some type. Non-zero of the eltype
               of the matrix is used to determine if a link exists,
               but currently no sense of magnitude
    C          Constant to fiddle with density of resulting layout
    MAXITER    Number of iterations we apply the forces
    INITTEMP   Initial "temperature", controls movement per iteration
"""
function layout_spring_adj{T}(adj_matrix::Array{T,2}; C=2.0, MAXITER=100, INITTEMP=2.0)

    size(adj_matrix, 1) != size(adj_matrix, 2) && error("Adj. matrix must be square.")
    const N = size(adj_matrix, 1)

    # Initial layout is random on the square [-1,+1]^2
    locs_x = 2*rand(N) .- 1.0
    locs_y = 2*rand(N) .- 1.0

    # The optimal distance bewteen vertices
    const K = C * sqrt(4.0 / N)

    # Store forces and apply at end of iteration all at once
    force_x = zeros(N)
    force_y = zeros(N)

    # Iterate MAXITER times
    @inbounds for iter = 1:MAXITER
        # Calculate forces
        for i = 1:N
            force_vec_x = 0.0
            force_vec_y = 0.0
            for j = 1:N
                i == j && continue
                d_x = locs_x[j] - locs_x[i]
                d_y = locs_y[j] - locs_y[i]
                d   = sqrt(d_x^2 + d_y^2)
                if adj_matrix[i,j] != zero(eltype(adj_matrix)) || adj_matrix[j,i] != zero(eltype(adj_matrix))
                    # F = d^2 / K - K^2 / d
                    F_d = d / K - K^2 / d^2
                else
                    # Just repulsive
                    # F = -K^2 / d^
                    F_d = -K^2 / d^2
                end
                # d  /          sin θ = d_y/d = fy/F
                # F /| dy fy    -> fy = F*d_y/d
                #  / |          cos θ = d_x/d = fx/F
                # /---          -> fx = F*d_x/d
                # dx fx
                force_vec_x += F_d*d_x
                force_vec_y += F_d*d_y
            end
            force_x[i] = force_vec_x
            force_y[i] = force_vec_y
        end
        # Cool down
        TEMP = INITTEMP / iter
        # Now apply them, but limit to temperature
        for i = 1:N
            force_mag  = sqrt(force_x[i]^2 + force_y[i]^2)
            scale      = min(force_mag, TEMP)/force_mag
            locs_x[i] += force_x[i] * scale
            #locs_x[i]  = max(-1.0, min(locs_x[i], +1.0))
            locs_y[i] += force_y[i] * scale
            #locs_y[i]  = max(-1.0, min(locs_y[i], +1.0))
        end
    end

    # Scale to unit square
    min_x, max_x = minimum(locs_x), maximum(locs_x)
    min_y, max_y = minimum(locs_y), maximum(locs_y)
    function scaler(z, a, b)
        2.0*((z - a)/(b - a)) - 1.0
    end
    map!(z -> scaler(z, min_x, max_x), locs_x, locs_x)
    map!(z -> scaler(z, min_y, max_y), locs_y, locs_y)

    return locs_x,locs_y
end
