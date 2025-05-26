using Yao, Yao.EasyBuild
using LinearAlgebra

function rabi_frequency(Ω_f::Float64, τ::Float64, α::Float64)
    T = τ/8
    a = exp(-(2*T)^2/(2*(α*T)^2))
    b = exp(-(τ/2-6*T)^2/(2*(α*T)^2))
    function O(t::Float64)
        if 0 < t <= 4*T
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

function calculate_pulse_area(Ω_f::Float64, τ::Float64, α::Float64)
    # 创建Rabi频率函数
    Ω(t) = Rabi_frequency(Ω_f, τ, α)(t)
    
    # 使用数值积分计算脉冲面积
    # 积分区间从 -8 到 τ
    dt = 0.001  # 积分步长
    t_range = 0:dt:τ
    area = 0.0
    
    for t in t_range
        area += abs(Ω(t)) * dt
    end
    
    return area
end

function effectiveHamiltonian1(Ω_0::Float64, Ω_1::Float64)

    ham = [0 0 Ω_0 Ω_0; 0 0 Ω_1 Ω_1; Ω_0 Ω_1 0 0; Ω_0 Ω_1 0 0]

    return ham
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
function initialize_superposition_state1()
    reg = zero_state(1, nlevel=4)
    H4 = create_four_level_hadamard()
    apply!(reg, put(1=>GeneralMatrixBlock(H4, nlevel=4)))
    #state_vector = zeros(ComplexF64, 64)
    #state_vector[1] = 1/2
    #state_vector[5] = 1/2
    #state_vector[17] = 1/2
    #state_vector[21] = 1/2
    #reg = ArrayReg(state_vector, nlevel=4)
    return reg
end

function evolve_eff1(reg, Ω_f::Float64, τ::Float64, α::Float64, Nt = 10000)
    Ω_0(t) = rabi_frequency(Ω_f, τ, α)(t)
    #/sqrt(2)
    Ω_1(t) = -1 * rabi_frequency(Ω_f, τ, α)(t)
    #/sqrt(2)
    state0 = Float64[]
    state1 = Float64[]
    stated = Float64[]
    statep = Float64[]
    dt = τ/Nt
    for i in 1:Nt
        t = (i-0.5)*dt
        ham_matrix = effectiveHamiltonian1(Ω_0(t), Ω_1(t))
        ham_block = GeneralMatrixBlock(ham_matrix;nlevel = 4)
        time_evo_op = time_evolve(ham_block,dt)
        apply!(reg, time_evo_op)
        state_0 = abs2(reg.state[1])
        state_1 = abs2(reg.state[2])
        state_d = abs2(reg.state[3])
        state_p = abs2(reg.state[4])
        push!(state0, state_0)
        push!(state1, state_1)
        push!(stated, state_d)
        push!(statep, state_p)
    end
    return state0, state1, stated, statep
end
