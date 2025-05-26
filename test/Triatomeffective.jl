using Test
using Plots
using LinearAlgebra
using Yao, Yao.EasyBuild
using ParityGate
import Yao: expect

@testset "ARP_pulse" begin
    Ω_f = 0.5069 * 2π
    Ω_c = 3.533 * 2π
    α = 0.5253
    τ = 3.0
    
    t_total = τ
    times = 0:0.01:τ
    t = times/τ
    # Calculate values for both functions
    Omega_0 = Rabi_frequency(Ω_f, τ, α)
    Omegavec = [Omega_0(t) for t in times]
    
    # Create the plot
    p = Plots.plot(t, [real(Omegavec)/2/pi], 
        label=["Rabi Frequency"],
        title="Gaussian Pulse of gates PE-X",
        xlabel="Time",
        ylabel="Amplitude",
        linewidth=2)
    display(p)
    @test p isa Plots.Plot
end

@testset "effectiveHamiltonian" begin
    Ω_0 = 0.5069 * 2π
    Ω_1 = 0.5069 * 2π


    ham_matrix = effectiveHamiltonian1(Ω_0, Ω_1)
    @test ham_matrix isa Matrix{Float64}
    
end

@testset "create_four_level_hadamard" begin
    hadamard = create_four_level_hadamard()
    @test hadamard isa Matrix{ComplexF64}
    @test size(hadamard) == (4, 4)
end


@testset "evolve_eff" begin
    Ω_f = 0.5069 * 2π
    Ω_c = 3.533 * 2π
    α = 0.529
    τ = 3.0
    J_12 = 0.01 * 2π
    J_dd = 0.005 * 2π
    Δ = 0.0

    ham_matrix = effectiveHamiltonian1(Ω_0, Ω_1)
    reg = initialize_superposition_state1()
    state0, state1, stated, statep = evolve_eff1(reg, Ω_f, τ, α)
    @test state0 isa Vector{Float64}
    @test state1 isa Vector{Float64}
    @test stated isa Vector{Float64}
    @test statep isa Vector{Float64}
end

@testset "paritygate" begin
    Ω_f = 0.5069 * 2π/2
    Ω_c = 3.533 * 2π
    α = 0.5253
    τ = 3.0
    Δ = 20 * Ω_c
    J_dd = Δ
    J_12 = 0.0
    Nt = 10000
    reg = zero_state(1, nlevel = 4)
    #reg = initialize_superposition_state1()
    print_table(reg)
    state0, state1, stated, statep = evolve_eff1(reg, Ω_f, τ, α)
    times = 0:τ/Nt:τ-τ/Nt
    p1 = Plots.plot(times, [state0, state1, stated, statep],
        label=["0" "1" "d" "p"],
        title="Parity gate",
        xlabel="Time",
        ylabel="Probability",
        linewidth=2)
    display(p1)
    @test p1 isa Plots.Plot
end