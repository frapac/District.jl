
function generate_db()
    println("Load data")
    db = @time load("$PATH_DATA/db.jld")
    getdb() = db
    return getdb
end

getdb = generate_db()

