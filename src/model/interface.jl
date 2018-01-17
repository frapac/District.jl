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
    price::Array{Float64}
    linker::AbstractConnection
end

swap!(p::PriceInterface, v::Array{Float64}) = copy!(p.price, v)


# Inteface for primal decomposition
struct FlowInterface <: AbstractInterface
    flow::Array{Float64}
end

swap!(p::FlowInterface, v::Array{Float64}) = p.flow = v
