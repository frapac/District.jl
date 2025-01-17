push!(LOAD_PATH, "..")


using District, StochDynamicProgramming

srand(10)

ALGO = "SDDP"

# load time span
ts = TimeSpan(150, 1)

################################################################################
# BUILDING MODEL
# initiate building
prof = District.ElecHouse(pv=20, heat=6, bat="bat0")
house = District.load(ts, prof)

# initial position
x0 = [.55, 2., 18., 18.]
# build SP model
District.build!(house, x0)


################################################################################
# SIMULATION
# Build simulator
nscen = 1000
sim = Simulator(house, nscen)

# Build policy
if ALGO == "MPC"
    pol = District.MPCPolicy(District.genforecast(ts, house.noises))
elseif ALGO == "SDDP"
    # Compute cuts :
    sddp = @time District.solve(house, SDDP(200))
    pol = District.HereAndNowDP(sddp.bellmanfunctions)
end

res = District.simulate(sim, pol)
