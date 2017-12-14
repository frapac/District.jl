
abstract type AbstractBuilding end


struct House
    name::Symbol
    devices::Vector{AbstractDevice}
    model
end

House()=House(:name, AbstractDevice[], nothing)



