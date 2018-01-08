push!(LOAD_PATH, "..")


using District, StochDynamicProgramming

ALGO = "MPC"

# Build problem
ts = TimeSpan(200, 1)
house = House(ts)

# Add devices
chp = MicroCHP("chp0")
add!(house, chp)

b = Battery("bat0")
add!(house, b)

hwt = ThermalHotWaterTank("twht0")
add!(house, hwt)

thm = R6C2("rt1988")
add!(house, thm)

heat = ThermalHeater(6.)
add!(house, heat)

District.link!(house, thm, heat)
District.link!(house, hwt, heat)

# Add noises
wdem = Demands(10, 1)
add!(house, wdem)


dynam = District.builddynamic(house)
load = District.buildload(house)

x0 = [.55, 2., 16., 16.]
District.build!(house, x0)


# Build simulator
sim = Simulator(house, 10)

# Build policy
if ALGO == "MPC"
    pol = District.MPCPolicy(District.genforecast(ts, house.noises))
elseif ALGO == "SDDP"
    # Compute cuts :
    sddp = @time District.solve(house, SDDP(100))
    pol = District.HereAndNowDP(sddp.bellmanfunctions)
end

res = District.simulate(sim, pol)
