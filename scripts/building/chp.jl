push!(LOAD_PATH, "..")


using District, StochDynamicProgramming

ALGO = "SDDP"

# Build problem
ts = TimeSpan(16, 1)
house = District.load(ts, District.CHPHouse(bat="bat0", nbins=10))

x0 = [.55, 5., 17., 17.]
x0 = [ 0.520956, 0, 18.0248, 17.5386]

District.build!(house, x0, info=:HD)

# Build simulator
sim = Simulator(house, 1000)
sim.model.controlCat[1] = :Cont

# Build policy
if ALGO == "MPC"
    pol = District.MPCPolicy(District.genforecast(ts, house.noises))
elseif ALGO == "SDDP"
    # Compute cuts:
    sddp = @time District.solve(house, SDDP(100))
    pol = District.HereAndNowDP(sddp.bellmanfunctions)
end

#= c, x, u = simulate(sddp, 100) =#
#= x0 = vec(mean(x[end, :, :], 1)) =#


res = District.simulate(sim, pol)
