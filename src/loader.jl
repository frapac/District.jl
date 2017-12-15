
function generate_db()
    println("Load data")
    db = @time load("$PATH_DATA/db.jld")
    return db
end

#= getdb = generate_db() =#


function loadweather()
    db = load("data/weather.jld")
    return db["weather"]["data"]
end
