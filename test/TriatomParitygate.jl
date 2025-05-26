using Test
using Plots
using LinearAlgebra
using Yao, Yao.EasyBuild
using ParityGate
import Yao: expect
using Symbolics

@testset "ARP_pulse" begin
    Ω_f = 0.5069 * 2π
    Ω_c = 3.533 * 2π
    α = 0.5253
    τ = 3.0
    
    t_total = τ
    times = 0:0.01:τ
    t = times/τ
    # Calculate values for both functions
    Omega_0 = rabi_frequency(Ω_f, τ, α)
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

    area = calculate_pulse_area(Ω_f, τ, α)
    println("pulse area: ", area)
end

@testset "effectiveHamiltonian" begin
    Ω_0 = 0.5069 * 2π
    Ω_1 = 0.5069 * 2π
    Ω_c = 3.533 * 2π
    α = 0.529
    τ = 3.0
    J_12 = 0.01 * 2π
    J_dd = 0.005 * 2π
    Δ = 0.0


    ham_matrix = effectiveHamiltonian(Ω_0, Ω_1, Ω_c, Δ, J_12, J_dd)
    @test ham_matrix isa Matrix{Float64}
    
end

@testset "create_four_level_hadamard" begin
    hadamard = create_four_level_hadamard()
    reg = initialize_superposition_state()
    print_table(reg)
    @test hadamard isa Matrix{ComplexF64}
    @test size(hadamard) == (4, 4)
end



@testset "paritygate" begin
    Ω_f = 0.5069 * 2π/2
    Ω_c = 3.533 * 2π/2
    α = 0.5253
    τ = 3.0
    Δ = 20 * Ω_c
    J_dd = Δ
    J_12 = 0.0
    Nt = 10000
    V_d = 66.978 * 2π
    V_pp = 3.7223 * 2π
    V_p = 238.23 * 2π

    reg = initialize_superposition_state()

    #print_table(reg)
    state000, state001, state110, state111, state010, state011, state100, state101 = evolve_eff(reg, Ω_f, τ, α, Ω_c, Δ, J_12, J_dd, V_p, V_pp, V_d, Nt)
    times = 0:τ/Nt:τ-τ/Nt
    p1 = Plots.plot(times, [state000, state001, state110, state111],
        label=["000" "001" "110" "111"],
        title="Parity gate",
        xlabel="Time",
        ylabel="Probability",
        linewidth=2)
    display(p1)
    p2 = Plots.plot(times, [state010, state011, state100, state101],
        label=["010" "011" "100" "101"],
        title="Parity gate",
        xlabel="Time",
        ylabel="Probability",
        linewidth=2)
    display(p2)
    #p3 = Plots.plot(times, [state000, state100, state011, state111],
     #   label=["000" "100" "011" "111"],
      #  title="Parity gate",
       # xlabel="Time",
       # ylabel="Probability",
       # linewidth=2)
    #   display(p3)
    #p4 = Plots.plot(times, [state110, state010, state101, state001],
     #   label=["110" "010" "101" "001"],
      #  title="Parity gate",
       # xlabel="Time",
       # ylabel="Probability",
       # linewidth=2)
    #  display(p4)
    @test p1 isa Plots.Plot
    @test p2 isa Plots.Plot
end

@testset "symbolic_test" begin
    ham = symbolic_hamiltonian()
# 打印矩阵
println("哈密顿量矩阵的符号形式：")

for i in 1:64
    for j in 1:64
        if !iszero(ham[i,j])
            println("H[$i,$j] = ", ham[i,j])
        end
    end
end

end
