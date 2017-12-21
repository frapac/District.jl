
################################################################################
# Decomposition schemes
################################################################################
@compat abstract type Decomposition end
immutable PriceDecomposition <: Decomposition end
immutable QuantDecomposition <: Decomposition end
immutable ADMMDecomposition <: Decomposition end



# Define atom for Bellman's update and simulation
@compat abstract type Atom end

struct PriceAtom <: Atom
    λ::Array
end

struct QuantAtom <: Atom
    F::Array
end

struct ADMMAtom <: Atom
    λ::Array
    ΔQ::Array
    penalty::Float64
end


################################################################################
# PROBLEM'S UPDATE
################################################################################

"""Update `JuMP.Model` objective with new multiplier `λ`."""
function updpb!(m::JuMP.Model, atom::PriceAtom, t, ny)
    λ = atom.λ[t, ny]
    pobj = m.ext[:obj]
    # we get the control, knowing that its last component is the allocation
    u = m[:u]
    # update objective
    @objective(m, :Min, pobj + λ*u[end])
end


"""Update `JuMP.Model` constraint with new allocation `f`."""
function updpb!(m::JuMP.Model, atom::QuantAtom, t, ny)
    f = atom.F[t, ny]
    # here, we assume that the flow f is unidimensional
    # we set the allocation in the JuMP Model
    JuMP.setRHS(m.ext[:alloc], f)
end

"""Update `JuMP.Model` constraint with augmented Lagrangian."""
function updpb!(m::JuMP.Model, atom::ADMMAtom, t, ny)
    τ = atom.penalty
    dq = atom.ΔQ[t, ny]
    λ = atom.λ[t, ny]
    pobj = m.ext[:obj]
    u = m[:u]
    # we define the augmented Lagrangian
    @objective(m, :Min, pobj + λ*(dq+u[end]) + τ*(dq + u[end])^2)
end
