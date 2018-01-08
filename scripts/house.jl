push!(LOAD_PATH, "..")


using District, StochDynamicProgramming

ALGO = "MPC"

# Build problem
ts = TimeSpan(200, 3)
house = House(ts)

# Add devices
b = Battery("bat0")
add!(house, b)

hwt = HotWaterTank("ehwt0")
add!(house, hwt)

thm = R6C2("rt1988")
add!(house, thm)


# Add noises
wdem = Demands(10, 1)
add!(house, wdem)
wpv = PVProduction(1, .15, 20, 0)
add!(house, wpv)

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
