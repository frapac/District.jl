# Generic functions to load data

# Load demands database
function getdb()
    db = load("data/db.jld")
    return db
end

# Load weather database
function loadweather()
    db = load("data/weather.jld")
    return db["weather"]["data"]
end
