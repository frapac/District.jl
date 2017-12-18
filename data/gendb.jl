# Read all CSV data and save them in HDF5 file

using JLD, HDF5

const PATH_DATA = "scenarios"
const PATH_HOUSE = "buildings"

if ARGS[1] == "weather"
    weather = readcsv("weather/weatherdata.csv", Float64, header=true)[1]

    jldopen("weather.jld", "w") do file
        g = g_create(file, "weather")
        g["data"] = weather
    end
else
    elec = readcsv("$PATH_DATA/elec200_15mn.csv", Float64)/1000
    occ  = readcsv("$PATH_DATA/occ200_15mn.csv", Float64)
    chp = Dict("De"=>elec,
            "Occ"=>occ)


    jldopen("db.jld", "w") do file

        g = g_create(file, "chp")
        g["De"]  = elec
        g["Occ"] = occ

        for i in 1:3
            g = g_create(file, "house$i")
            elec = readcsv("$PATH_HOUSE/house$i/elec.csv")[2:end, 2:end]/1000
            elec = Array{Float64, 2}(elec)
            mdhw = readcsv("$PATH_HOUSE/house$i/DHW.csv")[2:end, 2:end]
            mdhw = Array{Float64, 2}(mdhw)
            g["De"]  = elec
            g["DHW"] = mdhw
        end
    end

end
