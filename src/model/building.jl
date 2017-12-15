
abstract type AbstractBuilding end


struct House <: AbstractBuilding
    name::Symbol
    time::AbstractTimeSpan
    devices::Vector{AbstractDevice}
    noises::Vector{AbstractUncertainty}
    model
end

House()=House(:name, AbstractDevice[], nothing)


function build!(house::House) end



### DEPRECATED
function build_costs2()
    priceElec, Tcons = get_reference()
    count = 0
    function cost_model(m, t, x, u, w)
        zel1 = @JuMP.variable(m, lowerbound=0)
        @constraint(m, zel1 >= priceElec[t] * (w[1] + u[4] - w[5] - u[2] + u[1] + u[3]))

        zth1 = @JuMP.variable(m, lowerbound=0)
        @constraint(m, zth1 >= -Params.P_TH*(x[4] - Tcons[t] + 1))

        return JuMP.AffExpr(zel1 + zth1)
    end
    return cost_model
end


function final_cost(model, m)
    alpha = m[:alpha]
    #= w = JuMP.getvariable(m, :w) =#
    x = m[:x]
    u = m[:u]
    xf = m[:xf]
    @variable(m, cost)
    z1 = @JuMP.variable(m, lowerbound=0)
    @JuMP.constraint(m, z1 >= 6 - xf[2])
    @JuMP.constraint(m, alpha == PENAL_TANK*z1)
end


function final_cost_dh(model, m)
    alpha = m[:alpha]
    # get number of random noises
    ns = model.noises[end-1].supportSize

    x = m[:x]
    u = m[:u]
    xf = m[:xf]
    @JuMP.variable(m, z1[1:ns] >= 0)
    @JuMP.constraint(m, z1[i=1:ns] .>= 6 - xf[2, i])
    for i in 1:ns
        @JuMP.constraint(m, alpha[i] == PENAL_TANK*z1[i])
    end
end


# TODO: update
"""Load demands scenarios in memory."""
function get_aleas(ntime; scen="optim", idday=1, ndays=1)
    db = getdb()
    if scen == "optim"
        nf, nt = 1, 1000
    else
        nf, nt = 1001, 2000
    end

    elec = db["house1"]["De"]
    dhw  = db["house1"]["DHW"]
    aleas = zeros(ndays*96, 1000, 2)
    # Initial time :
    ti = (idday-1)*96 + 1
    # Final time :
    te = (idday + ndays - 1)*96
    # Select corresponding data :
    aleas[:, :, 1] = elec[ti:te, nf:nt]
    aleas[:, :, 2] = dhw[ti:te, nf:nt]
    return aleas
end


function build!(house::House)
#= function init_problem(;nclusts=N_MEDOIDS, idday=-1, x0=[.6, 6, 16, 16], =#
#=                       switch=false, ndays=1, tol_pv=TOL_PV, info=:HD, nrad=N_RAD) =#
    # Instantiate model:
    maxdhw = 20.
    heater = 6.

    xbounds = [(Params.BMIN, Params.BMAX), (0, 12), (-50, 100), (-50, 100)]
    ubounds = [(0, Params.DELTA_B_MAX), (0, Params.DELTA_B_MAX), (0., heater), (0, maxdhw)]

    # Set day to idday, otherwise Params.DAY as default value
    day = (idday > 0)? idday : Params.DAY
    weather = import_data(1, ndays, start=(day-1)*96+1)

    Pint = weather.radiative_internal_gain
    Pext = weather.radiative_external_gain
    α = cov(Pint, Pext)/var(Pext)
    corrpv = cov(weather.prod_pv, Pext)/var(Pext)

    text_f, phi_f = getweather(weather, ndays=ndays)

    demands = get_aleas(96, idday=(day-1)%7 + 1, ndays=ndays)
    laws = build_aleas(96, demands, text_f, phi_f, weather.prod_pv, nclusts,
                       nrad=nrad, epsrad=tol_pv)

    function dynamic(t, x, u, w)
        return [Params.ALPHA_B * x[1] + Params.DT*(Params.ρc*u[1] - 1/Params.ρd*u[2]), # Battery
                x[2] + Params.DT*(u[4]*Params.ETA_R - w[2]), # Hot water tank
                x[3] + Params.DT*3600/Params.Cw*(Params.Giw*(x[4]-x[3]) +
                                          Params.Gwe*(w[3]-x[3]) +
                                          Params.Re/(Params.Re+Params.Rw)*w[4] +
                                          Params.Ri/(Params.Ri+Params.Rs)*α*w[4] +
                                          1000*.35*u[3]), # wall's temperature
                x[4] + Params.DT*3600/Params.Ci*(Params.Giw*(x[3]-x[4]) +
                                          Params.Gie*(w[3]-x[4]) +
                                          Params.Rs/(Params.Ri+Params.Rs)*α*w[4] +
                                          1000*.65*u[3])] # inner temperature
    end

    fcost = (info == :HD)? final_cost : final_cost_dh
    model = StochDynamicProgramming.LinearSPModel(96*ndays, ubounds,
                                                  x0, build_costs2(),
                                                  dynamic, laws,
                                                  info=info,
                                                  Vfinal=fcost)

    println("Size of perturbations: ", laws[1].supportSize)
    set_state_bounds(model, xbounds)

    return model
end





"""Get real cost for assessment."""
function get_real_cost(day)
    priceElec, Tcons = get_reference()

    function real_cost(t, x, u, w)
        flow  = (w[1] + u[4]- w[5] - u[2] + u[1] + u[3])
        pelec = priceElec[t]*max(0, flow)

        temp  = -x[4] + Tcons[t] - 1
        pconfort = Params.P_TH * max(0, temp)

        return pelec + pconfort
    end

    return real_cost
end
realfinalcost(xf) = PENAL_TANK*max(0, 6 - xf[2])
