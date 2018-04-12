################################################################################
# District.jl
################################################################################
# Define interface between node and the rest of the graph
################################################################################

export FlowInterface, PriceInterface, QuadInterface, ZoneInterface

abstract type AbstractInterface end


# If no decomposition scheme, we use NoneInterface
struct NoneInterface <: AbstractInterface end
NoneInterface(values, linker) = NoneInterface()


# Interface for price decomposition
struct PriceInterface <: AbstractInterface
    values::Array{Float64}
    linker::AbstractConnection
end

swap!(p::PriceInterface, v::Array{Float64}) = copy!(p.values, v)


# Inteface for primal decomposition
struct FlowInterface <: AbstractInterface
    values::Array{Float64}
    linker::AbstractConnection
end

swap!(p::FlowInterface, v::Array{Float64}) = copy!(p.values, v)

# Interface for augmented Lagrangian
struct QuadInterface <: AbstractInterface
    τ::Float64
    values::Array{Float64}
    penal::Array{Float64}
    linker::AbstractConnection
end

swap!(p::QuadInterface, v::Array{Float64}) = copy!(p.values, v)
flow!(p::QuadInterface, f::Array{Float64}) = copy!(p.penal, f)

# Interface for zone
struct ZoneInterface <: AbstractInterface
    values::Array{Float64}
    # Matrix indicating position of border nodes
    linker::Array{Int64}
end

swap!(p::ZoneInterface, v::Array{Float64}) = copy!(p.values, v)
