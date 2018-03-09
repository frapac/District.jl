include("stress.jl")

"Get adjacence matrix of node-arc incidence matrix `A`."
function getadjacence(A::Array{Float64, 2})
    nnodes = size(A, 1)
    L = A * A'
    for i in 1:nnodes
        L[i, i] = 0
    end
    return -L
end

"Get average flow stored in `u`."
function getflow(u, operation=mean)
    nu = sum(District.ncontrols.(pb.nodes))
    q = u[:, :, nu+1:end]
    return vec(operation(operation(q, 2), 1))
end

"Display topology of graph specified by node-arc incidence matrix `A`."
function plotgraph(A)
    srand(11)
    nnodes, narcs = size(A)
    adjmat = getadjacence(A)
    posx, posy = layout_spring_adj(adjmat)

    scatter(posx, posy, s=100, color="k")

    for i = 2:nnodes, j=1:i
        if adjmat[i, j] != 0.
            plot([posx[i], posx[j]], [posy[i], posy[j]], c="k")
        end
    end
    axis("off")

end

"Display flow `q` on graph specified by node-arc incidence matrix `A`."
function plotflow(A, q)
    srand(11)
    nnodes, narcs = size(A)
    adjmat = getadjacence(A)
    posx, posy = layout_spring_adj(adjmat)


    qmax = maximum(abs.(q))

    idarc = 1
    for edge in 1:narcs
        pos = find(x->(x!=0), A[:, edge])
        α = abs(q[edge] / qmax)
        i, j = pos[1], pos[2]
        plot([posx[i], posx[j]], [posy[i], posy[j]], c=(α, 0., 1- α), lw=4., zorder=1)
    end
    scatter(posx, posy, s=150, color="k", zorder=2)
    axis("off")

end
