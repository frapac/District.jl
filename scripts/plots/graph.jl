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









"Display topology of graph specified by node-arc incidence matrix `A`."
function plotgraph(A)
    srand(11)
    figure()
    
    nnodes, narcs = size(A)
    adjmat = District.getadjacence(A)
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
function plotflow(A, q; darrow=false, offsetx=.12, offsety=.05)
    srand(11)
    nnodes, narcs = size(A)
    adjmat = District.getadjacence(A)
    posx, posy = layout_spring_adj(adjmat)


    qmax = maximum(abs.(q))

    figure()
    idarc = 1
    for edge in 1:narcs
        pos = find(x->(x!=0), A[:, edge])

        # find in and outgoing nodes (just matter of convention)
        if A[pos[1], edge] > 0
            i, j = pos[1], pos[2]
        else
            i, j = pos[2], pos[1]
        end

        α = abs(q[edge] / qmax)
        plot([posx[i], posx[j]], [posy[i], posy[j]], c=(α, 0., 1- α), lw=10*α, zorder=1)

        if darrow
            # start position
            x0, y0 = (q[edge] > 0.)? (posx[i], posy[i]) : (posx[j], posy[j])
            # stop position
            x1, y1 = (q[edge] > 0.)? (posx[j], posy[j]) : (posx[i], posy[i])

            # get vector coordinates
            dx = x1 - x0
            dy = y1 - y0
            # get norm
            ndxy = sqrt(dx^2 + dy^2)
            # build unit vectors
            u = [dx, dy] / ndxy
            v = [-dy, dx] / ndxy

            xstart, ystart = [x0, y0] + offsetx * u + offsety * v
            xstop, ystop = [x1, y1] - offsetx * u + offsety * v

            arrow(xstart, ystart, xstop - xstart, ystop - ystart,
                  head_width=0.03, head_length=0.05, fc="k")
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
    adjmat = District.getadjacence(A)
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
