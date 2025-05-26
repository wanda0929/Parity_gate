module ParityGate

using LinearAlgebra

include("TriatomParitygate.jl")
include("Triatomeffective.jl")

# Export your public API here
export Rabi_frequency, effectiveHamiltonian, evolve_eff, create_four_level_hadamard, initialize_superposition_state
export effectiveHamiltonian1, initialize_superposition_state1, evolve_eff1
export calculate_pulse_area
export symbolic_hamiltonian
export state_to_decimal, decimal_to_state
end # module ParityGate
