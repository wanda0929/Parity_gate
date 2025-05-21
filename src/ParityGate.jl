module ParityGate

using LinearAlgebra

include("3atomParitygate.jl")
include("3atomeffective.jl")

# Export your public API here
export Rabi_frequency, effectiveHamiltonian, evolve_eff, create_four_level_hadamard, initialize_superposition_state
export effectiveHamiltonian1, initialize_superposition_state1, evolve_eff1

end # module ParityGate
