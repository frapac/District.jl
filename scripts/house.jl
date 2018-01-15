push!(LOAD_PATH, "..")


using District, StochDynamicProgramming

ALGO = "MPC"

# load time span
ts = TimeSpan(200, 3)

################################################################################
# BUILDING MODEL
# initiate building
prof = District.ElecHouse(pv=20, heat=6, bat="bat0")
house = District.load(ts, prof)

# initial position
x0 = [.55, 2., 16., 16.]
# build SP model
District.build!(house, x0)


################################################################################
# SIMULATION
# Build simulator
nscen = 10
sim = Simulator(house, nscen)

# Build policy
if ALGO == "MPC"
    pol = District.MPCPolicy(District.genforecast(ts, house.noises))
elseif ALGO == "SDDP"
    # Compute cuts :
    sddp = @time District.solve(house, SDDP(100))
    pol = District.HereAndNowDP(sddp.bellmanfunctions)
end

res = District.simulate(sim, pol)
