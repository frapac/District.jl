push!(LOAD_PATH, "..")


using District, StochDynamicProgramming

ALGO = "MPC"

# Build problem
ts = TimeSpan(10, 1)
house = District.load(ts, District.CHPHouse(bat="bat0"))

x0 = [.55, 5., 17., 17.]
District.build!(house, x0, info=:HD)


# Build simulator
sim = Simulator(house, 10)

# Build policy
if ALGO == "MPC"
    pol = District.MPCPolicy(District.genforecast(ts, house.noises))
elseif ALGO == "SDDP"
    # Compute cuts :
    sddp = @time District.solve(house, SDDP(50))
    pol = District.HereAndNowDP(sddp.bellmanfunctions)
end

res = District.simulate(sim, pol)
