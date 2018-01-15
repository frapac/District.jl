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


function updpb!(m::JuMP.Model, p::PriceInterface, t, ny)
    # TODO: add ny
    λ = p.price[t]
    pobj = m.obj
    # we get the control, knowing that its last component is the allocation
    u = m[:u]
    # update objective
    @objective(m, :Min, pobj + λ*u[end])
end


# Inteface for primal decomposition
struct FlowInterface <: AbstractInterface
    flow::Array{Float64}
end

swap!(p::FlowInterface, v::Array{Float64}) = p.flow = v
