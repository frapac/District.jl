#using LightGraphs, GraphPlot, Compose, Colors
include("stress.jl")

hasbattery(node::House) = District.hasdevice(node, Battery)
haspv(node::House) = District.hasnoise(node, PVProduction)

function getshape(pb::Grid)
    shapes = String[]
    for d in pb.nodes
        if hasbattery(d) && haspv(d)
            push!(shapes, "darkred")
        elseif haspv(d)
            push!(shapes, "darkgreen")
        else
            push!(shapes, "black")
        end
    end
    shapes
end


"Get adjacence matrix of node-arc incidence matrix `A`."
function getadjacence(A::Array{Float64, 2})
    nnodes = size(A, 1)
    L = A * A'
    for i in 1:nnodes
        L[i, i] = 0
    end
    return -L
end

"Get laplacian matrix of node-arc incidence matrix `A` with weight q."
function getlaplacian(A::Array{Float64, 2}, q::Array{Float64, 1})
    nnodes = size(A, 1)
    narcs  = size(A, 2)

    Q = zeros(narcs,narcs)
    for i in 1:narcs
        Q[i,i] = q[i]
    end

    L = A * Q * A'

    return L
end

"Get average flow stored in `u`."
function getflow(u, operation=mean)
    nu = sum(District.ncontrols.(pb.nodes))
    q = abs.(u[:, :, nu+1:end])
    return vec(operation(operation(q, 2), 1))
end


"Display topology of graph specified by node-arc incidence matrix `A`."
function plotgraph(A)
    srand(11)
    figure()
    
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

"Display time-scenario mean absolute flow `q` on graph specified by node-arc incidence matrix `A`."
function plotflow(A, q; darrow=false, offset=.05, alpha=.1)
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
        plot([posx[i], posx[j]], [posy[i], posy[j]], c=(α, 0., 1- α), lw=10*α, zorder=1)

        if darrow
            # start position
            x0, y0 = (q[edge] > 0.)? (posx[i], posy[i]) : (posx[j], posy[j])
            # stop position
            x1, y1 = (q[edge] > 0.)? (posx[j], posy[j]) : (posx[i], posy[i])

            x0 -= offset; y0 -= offset
            x1 -= offset; y1 -= offset


            xstart = x0 + alpha*(x1 - x0)
            ystart = y0 + alpha*(y1 - y0)
            dx = (1 - 2*alpha) * (x1 - x0)
            dy = (1 - 2*alpha) * (y1 - y0)
            arrow(xstart, ystart, dx, dy,
                  head_width=0.01, head_length=0.05, fc="k")

        end
    end
    
    scatter(posx, posy, s=150, color=getshape(pb), zorder=2)


    axis("off")

end

"Display zonal decomposition and flow `q` on graph specified by node-arc incidence matrix `A`."
function plotzone(A, membership; darrow=false, offset=.05, alpha=.1, ncluster=3)


    srand(11)
    figure()
    nnodes, narcs = size(A)
    adjmat = getadjacence(A)
    posx, posy = layout_spring_adj(adjmat)


    qmax = maximum(q)

    idarc = 1
    for edge in 1:narcs
        pos = find(x->(x!=0), A[:, edge])
        α = abs(q[edge] / qmax)
        i, j = pos[1], pos[2]
        plot([posx[i], posx[j]], [posy[i], posy[j]], c=(α, 0., 1- α), lw=10*α, zorder=1 )

        if darrow
            # start position
            x0, y0 = (q[edge] > 0.)? (posx[i], posy[i]) : (posx[j], posy[j])
            # stop position
            x1, y1 = (q[edge] > 0.)? (posx[j], posy[j]) : (posx[i], posy[i])

            x0 -= offset; y0 -= offset
            x1 -= offset; y1 -= offset


            xstart = x0 + alpha*(x1 - x0)
            ystart = y0 + alpha*(y1 - y0)
            dx = (1 - 2*alpha) * (x1 - x0)
            dy = (1 - 2*alpha) * (y1 - y0)
            arrow(xstart, ystart, dx, dy,
                  head_width=0.06, head_length=0.05, fc="k", overhang=1)

        end
    end


    scatter(posx, posy, s=150, c=membership, zorder=2)
    axis("off")
    title("cluster")

end
