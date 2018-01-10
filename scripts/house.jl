push!(LOAD_PATH, "..")


using District, StochDynamicProgramming

ALGO = "MPC"

# Build problem
ts = TimeSpan(200, 3)

################################################################################
# BUILDING MODEL
# initiate building
house = House(ts)

# Add devices
devices = [Battery, ElecHotWaterTank, R6C2, ElecHeater]
names = ["bat0", "ehwt0", "rt1988", 6.]
dev = []
for (Device, name) in zip(devices, names)
    d = Device(name)
    add!(house, d)
    push!(dev, d)
end
bat, hwt, thm, heat = dev

# link heater to thermal envelope
District.link!(house, thm, heat)

# Add noises
# import demands
wdem = Demands(10, 1)
add!(house, wdem)
# import pv production
wpv = PVProduction(1, .15, 20, 0)
add!(house, wpv)

# link hot water tank with hot water demand
District.link!(house, hwt, wdem)

# build objective: we penalize elec and thermal comfort.
set!(house, EDFPrice(ts))
set!(house, ComfortPrice(ts))



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
