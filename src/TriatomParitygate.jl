using Yao, Yao.EasyBuild
using LinearAlgebra
using Symbolics
function rabi_frequency(Ω_f::Float64, τ::Float64, α::Float64)
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

function state_to_decimal(state::String)
    # 将能级标记转换为数字
    level_to_num = Dict('0'=>0, '1'=>1, 'd'=>2, 'p'=>3)
    
    # 确保输入是三位
    if length(state) != 3
        error("状态必须是三位")
    end
    
    # 转换每一位
    a = level_to_num[state[1]]  # 最高位（第一个原子）
    b = level_to_num[state[2]]  # 中间位（第二个原子）
    c = level_to_num[state[3]]  # 最低位（第三个原子）
    
    # 转换为十进制
    # 注意：这里要加1因为Yao中的索引从1开始
    decimal = 1 + a*16 + b*4 + c
    
    return decimal
end

# 反向转换：从十进制到态
function decimal_to_state(decimal::Int)
    # 减1因为Yao中的索引从1开始
    n = decimal - 1
    
    # 转换为四进制
    a = div(n, 16)  # 最高位
    b = div(mod(n, 16), 4)  # 中间位
    c = mod(n, 4)  # 最低位
    
    # 将数字转换为能级标记
    num_to_level = Dict(0=>'0', 1=>'1', 2=>'d', 3=>'p')
    
    return string(num_to_level[a], num_to_level[b], num_to_level[c])
end

function effectiveHamiltonian(Ω_0::Float64, Ω_1::Float64, Ω_c::Float64, Δ::Float64, J_12::Float64, J_dd::Float64, V_p::Float64, V_pp::Float64, V_d::Float64)
    #ham_1c = [0 0 Ω_c 0; 0 0 0 Ω_c; Ω_c 0 -Δ  0; 0 Ω_c 0 Δ]
    ham_1c = [-Δ/2 0 0 Ω_c; 0 Δ/2 Ω_c 0; 0 Ω_c -Δ/2 0; Ω_c 0 0 Δ/2]
    ham_3t = [0 0 Ω_0 Ω_0; 0 0 Ω_1 Ω_1; Ω_0 Ω_1 0 0; Ω_0 Ω_1 0 0]
    identity = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]
    
    # 初始化总哈密顿量
#    ham_matrix = zeros(ComplexF64, 64, 64)  # 假设是 6 量子比特系统
 #   ham_tot = kron(ham_1c, identity, identity) + kron(identity, ham_1c, identity) + kron(identity, identity, ham_3t)
  #  ham_matrix[1:64, 1:64] = ham_tot


    ham_matrix = kron(ham_1c, identity, identity) + kron(identity, ham_1c, identity) + kron(identity, identity, ham_3t)
    
    # 在矩阵添加元素
    ham_matrix[state_to_decimal("ddp"), state_to_decimal("pdd")] = J_dd
    ham_matrix[state_to_decimal("pdd"), state_to_decimal("ddp")] = J_dd
    ham_matrix[state_to_decimal("dpp"), state_to_decimal("ppd")] = J_dd
    ham_matrix[state_to_decimal("ppd"), state_to_decimal("dpp")] = J_dd
    ham_matrix[state_to_decimal("ddp"), state_to_decimal("dpd")] = J_dd
    ham_matrix[state_to_decimal("dpd"), state_to_decimal("ddp")] = J_dd
    ham_matrix[state_to_decimal("pdp"), state_to_decimal("ppd")] = J_dd
    ham_matrix[state_to_decimal("ppd"), state_to_decimal("pdp")] = J_dd

    ham_matrix[state_to_decimal("dpd"), state_to_decimal("pdd")] = J_12
    ham_matrix[state_to_decimal("pdd"), state_to_decimal("dpd")] = J_12
    ham_matrix[state_to_decimal("dpp"), state_to_decimal("pdp")] = J_12
    ham_matrix[state_to_decimal("pdp"), state_to_decimal("dpp")] = J_12

    ham_matrix[state_to_decimal("dp1"), state_to_decimal("pd1")] = J_12
    ham_matrix[state_to_decimal("pd1"), state_to_decimal("dp1")] = J_12
    ham_matrix[state_to_decimal("dp0"), state_to_decimal("pd0")] = J_12
    ham_matrix[state_to_decimal("pd0"), state_to_decimal("dp0")] = J_12
    
    ham_matrix[state_to_decimal("d1p"), state_to_decimal("p1d")] = J_dd
    ham_matrix[state_to_decimal("p1d"), state_to_decimal("d1p")] = J_dd
    ham_matrix[state_to_decimal("d0p"), state_to_decimal("p0d")] = J_dd
    ham_matrix[state_to_decimal("p0d"), state_to_decimal("d0p")] = J_dd

    ham_matrix[state_to_decimal("1dp"), state_to_decimal("1pd")] = J_dd
    ham_matrix[state_to_decimal("1pd"), state_to_decimal("1dp")] = J_dd
    ham_matrix[state_to_decimal("0dp"), state_to_decimal("0pd")] = J_dd
    ham_matrix[state_to_decimal("0pd"), state_to_decimal("0dp")] = J_dd

    ham_matrix[state_to_decimal("p0p"), state_to_decimal("p0p")] = V_p
    ham_matrix[state_to_decimal("p1p"), state_to_decimal("p1p")] = V_p
    ham_matrix[state_to_decimal("pdp"), state_to_decimal("pdp")] = V_p
    ham_matrix[state_to_decimal("ppp"), state_to_decimal("ppp")] = 2 * V_p + V_pp

    ham_matrix[state_to_decimal("0pp"), state_to_decimal("0pp")] = V_d
    ham_matrix[state_to_decimal("1pp"), state_to_decimal("1pp")] = V_d
    ham_matrix[state_to_decimal("dpp"), state_to_decimal("dpp")] = V_d

    ham_matrix[state_to_decimal("d0d"), state_to_decimal("d0d")] = V_d
    ham_matrix[state_to_decimal("d1d"), state_to_decimal("d1d")] = V_d
    ham_matrix[state_to_decimal("dpd"), state_to_decimal("dpd")] = V_d
    ham_matrix[state_to_decimal("ddd"), state_to_decimal("ddd")] = 2 * V_d

    ham_matrix[state_to_decimal("0dd"), state_to_decimal("0dd")] = V_d
    ham_matrix[state_to_decimal("1dd"), state_to_decimal("1dd")] = V_d
    ham_matrix[state_to_decimal("pdd"), state_to_decimal("pdd")] = V_d




    

    #ham_matrix[12, 15] = J_dd
    #ham_matrix[15, 12] = J_dd

    #ham_matrix[28, 31] = J_dd
    #ham_matrix[31, 28] = J_dd

    #ham_matrix[36, 51] = J_dd
    #ham_matrix[51, 36] = J_dd

    #ham_matrix[40, 55] = J_dd
    #ham_matrix[55, 40] = J_dd

    #ham_matrix[45, 57] = J_12
    #ham_matrix[57, 45] = J_12

    #ham_matrix[46, 58] = J_12
    #ham_matrix[58, 46] = J_12

    #ham_matrix[44, 47] = J_dd
    #ham_matrix[47, 44] = J_dd

    #ham_matrix[44, 59] = J_dd
    #ham_matrix[59, 44] = J_dd

    #ham_matrix[48, 63] = J_dd
    #ham_matrix[63, 48] = J_dd

    #ham_matrix[60, 63] = J_dd
    #ham_matrix[63, 60] = J_dd

    #ham_matrix[59, 47] = J_12
    #ham_matrix[47, 59] = J_12

    #ham_matrix[48, 60] = J_12
    #ham_matrix[60, 48] = J_12

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

function evolve_eff(reg, Ω_f::Float64, τ::Float64, α::Float64, Ω_c::Float64, Δ::Float64, J_12::Float64, J_dd::Float64, V_p::Float64, V_pp::Float64, V_d::Float64, Nt = 10000)
    state000 = Float64[]
    state001 = Float64[]
    state110 = Float64[]
    state111 = Float64[]
    state010 = Float64[]
    state011 = Float64[]
    state100 = Float64[]
    state101 = Float64[]
    Ω_0(t) = rabi_frequency(Ω_f, τ, α)(t)
    #/sqrt(2)
    Ω_1(t) = - rabi_frequency(Ω_f, τ, α)(t)
    #/sqrt(2)
    dt = τ/Nt
    for i in 1:Nt
        t = (i-0.5)*dt
        ham_matrix = effectiveHamiltonian(Ω_0(t), Ω_1(t), Ω_c, Δ, J_12, J_dd, V_p, V_pp, V_d)
        ham_block = GeneralMatrixBlock(ham_matrix;nlevel = 4)
        time_evo_op = time_evolve(ham_block,dt)
        apply!(reg, time_evo_op)
        state_000 = abs2(reg.state[state_to_decimal("000")])
        state_001 = abs2(reg.state[state_to_decimal("001")])
        state_110 = abs2(reg.state[state_to_decimal("110")])
        state_111 = abs2(reg.state[state_to_decimal("111")])
        state_010 = abs2(reg.state[state_to_decimal("010")])
        state_011 = abs2(reg.state[state_to_decimal("011")])
        state_100 = abs2(reg.state[state_to_decimal("100")])
        state_101 = abs2(reg.state[state_to_decimal("101")])
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

function symbolic_hamiltonian()
    # 定义符号变量
    @variables Ω_0 Ω_1 Ω_c Δ J_12 J_dd
    
    # 创建基本矩阵
    
    ham_1c = [-Δ/2 0 0 Ω_c; 0 Δ/2 Ω_c 0; 0 Ω_c -Δ/2 0; Ω_c 0 0 Δ/2]
    #ham_1c = [0 0 0 Ω_c; 0 0 Ω_c 0; 0 Ω_c Δ 0; Ω_c 0 0 -Δ]
    ham_3t = [0 0 Ω_0 Ω_0; 0 0 Ω_1 Ω_1; Ω_0 Ω_1 0 0; Ω_0 Ω_1 0 0]
    identity = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]
    
    # 计算总哈密顿量
    ham_tot = kron(ham_1c, identity, identity) + kron(identity, ham_1c, identity) + kron(identity, identity, ham_3t)
    
    # 添加耦合项
    ham_matrix = copy(ham_tot)
    
    # 在矩阵添加元素
    ham_matrix[state_to_decimal("ddp"), state_to_decimal("pdd")] = J_dd
    ham_matrix[state_to_decimal("pdd"), state_to_decimal("ddp")] = J_dd
    ham_matrix[state_to_decimal("dpp"), state_to_decimal("ppd")] = J_dd
    ham_matrix[state_to_decimal("ppd"), state_to_decimal("dpp")] = J_dd
    ham_matrix[state_to_decimal("ddp"), state_to_decimal("dpd")] = J_dd
    ham_matrix[state_to_decimal("dpd"), state_to_decimal("ddp")] = J_dd
    ham_matrix[state_to_decimal("pdp"), state_to_decimal("ppd")] = J_dd
    ham_matrix[state_to_decimal("ppd"), state_to_decimal("pdp")] = J_dd

    ham_matrix[state_to_decimal("dpd"), state_to_decimal("pdd")] = J_12
    ham_matrix[state_to_decimal("pdd"), state_to_decimal("dpd")] = J_12
    ham_matrix[state_to_decimal("dpp"), state_to_decimal("pdp")] = J_12
    ham_matrix[state_to_decimal("pdp"), state_to_decimal("dpp")] = J_12

    ham_matrix[state_to_decimal("dp1"), state_to_decimal("pd1")] = J_12
    ham_matrix[state_to_decimal("pd1"), state_to_decimal("dp1")] = J_12
    ham_matrix[state_to_decimal("dp0"), state_to_decimal("pd0")] = J_12
    ham_matrix[state_to_decimal("pd0"), state_to_decimal("dp0")] = J_12
    
    ham_matrix[state_to_decimal("d1p"), state_to_decimal("p1d")] = J_dd
    ham_matrix[state_to_decimal("p1d"), state_to_decimal("d1p")] = J_dd
    ham_matrix[state_to_decimal("d0p"), state_to_decimal("p0d")] = J_dd
    ham_matrix[state_to_decimal("p0d"), state_to_decimal("d0p")] = J_dd

    ham_matrix[state_to_decimal("1dp"), state_to_decimal("1pd")] = J_dd
    ham_matrix[state_to_decimal("1pd"), state_to_decimal("1dp")] = J_dd
    ham_matrix[state_to_decimal("0dp"), state_to_decimal("0pd")] = J_dd
    ham_matrix[state_to_decimal("0pd"), state_to_decimal("0dp")] = J_dd
        
    return ham_matrix
end

