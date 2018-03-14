################################################################################
# District.jl
################################################################################
# Definition of our toy graph to test zonal DADP
################################################################################

# we build twelve houses

houseArray = Vector{House}(nnodes)
houseArray[1] = load(ts, ElecHouse(pv=10, bat="", nbins=nbins, idhouse=1))
houseArray[2] = load(ts, ElecHouse(pv=10, bat="", nbins=nbins, idhouse=2))
houseArray[3] = load(ts, ElecHouse(pv=10, bat="", nbins=nbins, idhouse=3))

houseArray[4] = load(ts, ElecHouse(pv=0, bat="bat0", nbins=nbins, idhouse=3))
houseArray[5] = load(ts, ElecHouse(pv=0, bat="", nbins=nbins, idhouse=2))
houseArray[6] = load(ts, ElecHouse(pv=4, bat="", nbins=nbins, idhouse=1))

houseArray[7] = load(ts, ElecHouse(pv=0, bat="bat0", nbins=nbins, idhouse=1))
houseArray[8] = load(ts, ElecHouse(pv=0, bat="", nbins=nbins, idhouse=2))
houseArray[9] = load(ts, ElecHouse(pv=0, bat="", nbins=nbins, idhouse=3))

houseArray[10] = load(ts, ElecHouse(pv=0, bat="", nbins=nbins, idhouse=3))
houseArray[11] = load(ts, ElecHouse(pv=0, bat="", nbins=nbins, idhouse=1))
houseArray[12] = load(ts, ElecHouse(pv=4, bat="bat0", nbins=nbins, idhouse=2))



iniArray = Array{Array{Float64,1},1}(nnodes)
for i in 1:nnodes
    if size(houseArray[i].devices,1) >= 4
        iniArray[i] = [.55, 2., 20., 20.]
    else
        iniArray[i] = [2., 20., 20.]
    end
end

xini = Dict(houseArray[i]=> iniArray[i] for i in 1:nnodes)

# Define network
net = Network(ts, A)
net.k2 = 1e-2 # transport cost parameters
net.k1 = 1e-3
# Define corresponding grid
pb = District.Grid(ts, houseArray, net)