# import numpy as np
# import matplotlib.pyplot as plt
# from scipy.sparse import diags
# from scipy.integrate import solve_ivp

# # 参数设置
# L = 200              # 格点总数
# i0 = 100             # 初始粒子所在位置
# t_span = [0, 50]     # 时间范围
# t_eval = [1, 10, 20, 50]  # 要绘图的特定时间点

# def hamiltonian(L):
#     """构造哈密顿量 H (稀疏矩阵形式)"""
#     diagonals = [-np.ones(L), np.ones(L)]  # 主对角线和偏移对角线
#     H = diags(diagonals, [-1, 1], shape=(L, L), format='csr')
#     H[0, -1] = H[-1, 0] = 1  # 周期边界条件
#     return -H  # 哈密顿量

# # 初始条件
# def initial_wavefunction(L, i0):
#     psi0 = np.zeros(L, dtype=complex)
#     psi0[i0] = 1.0  # 粒子位于 i0 处
#     return psi0

# # 薛定谔方程右侧的微分函数
# def schrodinger_eq(t, psi, H):
#     return -1j * H @ psi

# # 数值求解
# def solve_schrodinger(L, i0, t_span, t_eval):
#     H = hamiltonian(L)
#     psi0 = initial_wavefunction(L, i0)
#     sol = solve_ivp(schrodinger_eq, t_span, psi0, t_eval=t_eval, args=(H,), method='RK45')
#     return sol.y  # 返回波函数

# # 绘图函数
# def plot_density(psi_values, t_eval):
#     plt.figure(figsize=(10, 6))
#     for idx, t in enumerate(t_eval):
#         density = np.abs(psi_values[:, idx])**2
#         plt.plot(density, label=f't = {t}')
#     plt.title('Particle Density Distribution at Different Times')
#     plt.xlabel('Lattice Site i')
#     plt.ylabel("hhh")
#     plt.legend()
#     plt.grid(True)
#     plt.show()

# # 主程序
# psi_values = solve_schrodinger(L, i0, t_span, t_eval)
# plot_density(psi_values, t_eval)


import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags
from scipy.integrate import solve_ivp

# 参数设置
L = 200           # 总格点数
i0 = 100          # 初始粒子位置
time_points = [1, 10, 20, 50]  # 观察的时间点
t_span = (0, 50)  # 时间区间

# 构造哈密顿量 (Hamiltonian)
def hamiltonian(L):
    diagonals = [-1 * np.ones(L), np.ones(L)]  # 对角线值
    H = diags(diagonals, offsets=[-1, 1], shape=(L, L), format='csr')
    # 周期边界条件
    H[0, -1] = H[-1, 0] = 1
    return -H

# 初始波函数
def initial_wavefunction(L, i0):
    psi_0 = np.zeros(L, dtype=complex)
    psi_0[i0] = 1.0  # 粒子初始位置
    return psi_0

# 薛定谔方程的右端项 dψ/dt = -i H ψ
def schrodinger_eq(t, psi, H):
    return -1j * H @ psi

# 数值求解
def solve_schrodinger(L, i0, t_span, time_points):
    H = hamiltonian(L)  # 哈密顿量
    psi_0 = initial_wavefunction(L, i0)  # 初始波函数
    solution = solve_ivp(schrodinger_eq, t_span, psi_0, t_eval=time_points, args=(H,))
    return np.abs(solution.y)**2  # 计算密度分布 |ψ|²

# 绘制每个时刻的粒子密度分布单独图表
def plot_single_time_density(density, time_points):
    for idx, t in enumerate(time_points):
        plt.figure(figsize=(8, 6))
        plt.plot(density[:, idx], label=f't = {t}')
        plt.title(f"粒子密度分布 |ψ|² 在 t = {t} 时刻")
        plt.xlabel("格点编号 i")
        plt.ylabel("密度分布 ρ_i(t)")
        plt.legend()
        plt.grid(True)
        plt.show()

# 计算并绘制每个时刻单独图表
density = solve_schrodinger(L, i0, t_span, time_points)
plot_single_time_density(density, time_points)
