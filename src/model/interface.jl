################################################################################
# District.jl
################################################################################
# Define interface between node and the rest of the graph
################################################################################

export FlowInterface, PriceInterface

abstract type AbstractInterface end


# If no decomposition scheme, we use NoneInterface
struct NoneInterface <: AbstractInterface end


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
