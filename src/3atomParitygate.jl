using Yao, Yao.EasyBuild
using LinearAlgebra

function Rabi_frequency(Ω_f::Float64, τ::Float64, α::Float64)
    T = τ/8
    a = exp(-(2*T)^2/(2*(α*T)^2))
    b = exp(-(τ/2-6*T)^2/(2*(α*T)^2))
    function O(t::Float64)
        if -8 < t <= 4*T
            return Ω_f * (exp(-(t-2*T)^2/(2*(α*T)^2))-a)/(1-a)
            #return real(Ω * (exp(-(t+5)^2/(2*τ^2)-im * α * (t+5)^2/2)))
        elseif 4*T < t <= τ
            return Ω_f * (exp(-(t-6*T)^2/(2*(α*T)^2))-b)/(1-b)
        else
            return 0.0
        end
    end
    return O
end

function effectiveHamiltonian(Ω_0::Float64, Ω_1::Float64, Ω_c::Float64, Δ::Float64, J_12::Float64, J_dd::Float64)
    ham_1c = [0 0 Ω_c 0; 0 0 0 Ω_c; Ω_c 0 -Δ  0; 0 Ω_c 0 Δ]
    ham_3t = [0 0 Ω_0 Ω_0; 0 0 Ω_1 Ω_1; Ω_0 Ω_1 0 0; Ω_0 Ω_1 0 0]
    identity = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]
    
    # 初始化总哈密顿量
#    ham_matrix = zeros(ComplexF64, 64, 64)  # 假设是 6 量子比特系统
 #   ham_tot = kron(ham_1c, identity, identity) + kron(identity, ham_1c, identity) + kron(identity, identity, ham_3t)
  #  ham_matrix[1:64, 1:64] = ham_tot


    ham_matrix = kron(ham_1c, identity, identity) + kron(identity, ham_1c, identity) + kron(identity, identity, ham_3t)
    
    # 在矩阵添加元素
    ham_matrix[12, 15] = J_dd
    ham_matrix[15, 12] = J_dd

    ham_matrix[28, 31] = J_dd
    ham_matrix[31, 28] = J_dd

    ham_matrix[36, 51] = J_dd
    ham_matrix[51, 36] = J_dd

    ham_matrix[40, 55] = J_dd
    ham_matrix[55, 40] = J_dd

    ham_matrix[45, 57] = J_12
    ham_matrix[57, 45] = J_12

    ham_matrix[46, 58] = J_12
    ham_matrix[58, 46] = J_12

    ham_matrix[44, 47] = J_dd
    ham_matrix[47, 44] = J_dd

    ham_matrix[44, 59] = J_dd
    ham_matrix[59, 44] = J_dd

    ham_matrix[48, 63] = J_dd
    ham_matrix[63, 48] = J_dd

    ham_matrix[60, 63] = J_dd
    ham_matrix[63, 60] = J_dd

    ham_matrix[59, 47] = J_12
    ham_matrix[47, 59] = J_12

    ham_matrix[48, 60] = J_12
    ham_matrix[60, 48] = J_12

    return ham_matrix
end

"""
    create_four_level_hadamard()

创建四能级系统的Hadamard门。这个门会将|0⟩态转换为(|0⟩ + |1⟩)/√2，
同时保持|r1⟩和|r2⟩态不变。

# Returns
一个4x4的矩阵，表示四能级系统的Hadamard门
"""
function create_four_level_hadamard()
    # 创建四能级Hadamard门
    H4 = zeros(ComplexF64, 4, 4)
    H4[1:2, 1:2] = [1 1; 1 -1] / sqrt(2)  # 对|0⟩和|1⟩应用标准Hadamard门
    H4[3, 3] = 1  # 保持|r1⟩不变
    H4[4, 4] = 1  # 保持|r2⟩不变
    return H4
end

"""
    initialize_superposition_state()

初始化量子态，其中第一和第二个原子处于 |0⟩ 和 |1⟩ 的叠加态，第三个原子处于 |0⟩ 态。
具体来说，第一和第二个原子的态为 (|0⟩ + |1⟩)/√2。

# Returns
一个包含指定叠加态的量子寄存器
"""
function initialize_superposition_state()
    reg = zero_state(3, nlevel=4)
    H4 = create_four_level_hadamard()
    apply!(reg, put(3=>GeneralMatrixBlock(H4, nlevel=4)))
    apply!(reg, put(2=>GeneralMatrixBlock(H4, nlevel=4)))

    #state_vector = zeros(ComplexF64, 64)
    #state_vector[1] = 1/2
    #state_vector[5] = 1/2
    #state_vector[17] = 1/2
    #state_vector[21] = 1/2
    #reg = ArrayReg(state_vector, nlevel=4)
    return reg
end

function evolve_eff(reg, Ω_f::Float64, τ::Float64, α::Float64, Ω_c::Float64, Δ::Float64, J_12::Float64, J_dd::Float64, Nt = 10000)
    state000 = Float64[]
    state001 = Float64[]
    state110 = Float64[]
    state111 = Float64[]
    state010 = Float64[]
    state011 = Float64[]
    state100 = Float64[]
    state101 = Float64[]
    Ω_0(t) = Rabi_frequency(Ω_f, τ, α)(t)/sqrt(2)
    Ω_1(t) = - Rabi_frequency(Ω_f, τ, α)(t)/sqrt(2)
    dt = τ/Nt
    for i in 1:Nt
        t = (i-0.5)*dt
        ham_matrix = effectiveHamiltonian(Ω_0(t), Ω_1(t), Ω_c, Δ, J_12, J_dd)
        ham_block = GeneralMatrixBlock(ham_matrix;nlevel = 4)
        time_evo_op = time_evolve(ham_block,dt)
        apply!(reg, time_evo_op)
        state_000 = abs2(reg.state[1])
        state_001 = abs2(reg.state[2])
        state_110 = abs2(reg.state[21])
        state_111 = abs2(reg.state[22])
        state_010 = abs2(reg.state[5])
        state_011 = abs2(reg.state[6])
        state_100 = abs2(reg.state[17])
        state_101 = abs2(reg.state[18])
        push!(state000, state_000)
        push!(state001, state_001)
        push!(state110, state_110)
        push!(state111, state_111)
        push!(state010, state_010)
        push!(state011, state_011)
        push!(state100, state_100)
        push!(state101, state_101)
    end
    return state000, state001, state110, state111, state010, state011, state100, state101
end
