
function getdb()
    db = load("data/db.jld")
    return db
end



function loadweather()
    db = load("data/weather.jld")
    return db["weather"]["data"]
end
