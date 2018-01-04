################################################################################
# District.jl
################################################################################
# Compute solar irradiation.
# - Currently, compute irradiation only for walls with fixed orientation.
# - Irradiation is computed with Direct and Diffuse horizontal irradiation
#   (BHI and DHI), with the Perez model.
################################################################################

#TODO: rename get_irradiation
export get_irradiation

"""Get declinaison."""
get_declinaison(day::Int) = asin(0.398 * sin((0.985 * day - 80)*pi/180))

"""Get solar time."""
get_time_theta(t::Float64) = 15*(t - 12) * pi/180


"""Get azimuth of sun for given time, day and altitude."""
function get_azimuth(t, day, altitude)
    declinaison = get_declinaison(day)
    horaire = get_time_theta(t)
    lambda = asin(cos(declinaison) * cos(horaire) / cos(altitude))
    return lambda
end


"""Get zenith angle for given time, day, and latitude."""
function get_zenith(t::Float64, day::Int, latitude)
    declinaison = get_declinaison(day)
    Ah = get_time_theta(t)
    return asin(sin(latitude) * sin(declinaison) + cos(latitude) * cos(declinaison) * cos(Ah))
end


""" Get solar irradiation."""
function get_irradiation(env::R6C2, ts::TimeSpan)
    day = ts.day
    ntime = ntimesteps(ts)
    DT = ts.δt

    irradbeam = loadweather(BHI(), ts)
    irraddiffuse = loadweather(DHI(), ts)

    # orientation of walls
    beta = [pi/2, pi/2, pi/2, pi/2, 0]
    gamma = [-pi, -pi/2, 0, pi/2, 0]

    G_extrad = zeros(ntime)
    G_intrad = zeros(ntime)

    for t=1:ntime
        Ib = irradbeam[t]
        Id = irraddiffuse[t]

        theta_zenith = get_zenith(t*DT, day, env.latitude)
        azimut = get_azimuth(t*DT, day, env.altitude)
        incidence = cos.(theta_zenith - beta).*cos.(azimut - gamma)

        Ib_vec = max.(0, Ib * incidence)
        Id_vec = Id .* (1 + cos.(beta))/2
        Ig = Ib .* cos(theta_zenith) + Id
        Ir_vec = env.albedo * Ig .* (1 - cos.(beta))/2

        Itot = Ib_vec + Id_vec + Ir_vec
        G_extrad[t] = sum(Itot.*env.Surf_wall)*env.esp
        G_intrad[t] = sum(Itot*env.fframe.*env.Surf_window)*env.Fv
    end

    return G_extrad, G_intrad
end
