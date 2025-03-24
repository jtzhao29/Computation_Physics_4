import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import root_scalar

# 物理常量
hbar = 1.0   # 普朗克常数
m = 1.0      # 质量
x_max = 10   # 空间范围
N = 500      # 离散步数

x = np.linspace(-x_max, x_max, N)  # 定义空间网格
dx = x[1] - x[0]                  # 网格步长

# 势能函数定义
def potential(x):
    return x**2 + 0.25 * x**4

V = potential(x)

# 龙格-库塔法积分薛定谔方程
def schrodinger(t, y, E):
    # y[0] = φ(t), y[1] = φ'(t)
    phi, dphi = y
    idx_t = int(min(max(t, 0), N - 1))  # 确保索引在有效范围内
    ddphi = -(2 * m / hbar**2) * (V[idx_t] - E) * phi
    return [dphi, ddphi]

# 匹配函数 f(E)
def matching_condition(E):
    # 初始条件
    y_left = [1e-3, 1e-2]  # 左侧初值 φ 和 φ'
    y_right = [1e-3, -1e-2]  # 右侧初值 φ 和 φ'
    
    # 左侧积分到 x_r
    sol_left = solve_ivp(schrodinger, [0, N // 2], y_left, args=(E,),
                         t_eval=np.linspace(0, N // 2, 100))
    # 右侧积分到 x_r
    sol_right = solve_ivp(schrodinger, [N - 1, N // 2], y_right, args=(E,),
                          t_eval=np.linspace(N - 1, N // 2, 100)[::-1])
    
    # 匹配点检查
    phi_left = sol_left.y[0, -1]
    dphi_left = sol_left.y[1, -1]
    phi_right = sol_right.y[0, -1]
    dphi_right = sol_right.y[1, -1]
    
    f = (phi_right - phi_left) - (dphi_right - dphi_left)
    return f

# # 调试匹配函数
# E_test = np.linspace(0.1, 20, 100)
# f_values = [matching_condition(E) for E in E_test]

# plt.plot(E_test, f_values)
# plt.axhline(0, color='red', ls='--')  # 添加零轴作为参考
# plt.xlabel("E")
# plt.ylabel("f(E)")
# plt.title("匹配函数 f(E) 的分布")
# plt.grid()
# plt.show()

# 使用动态调整区间查找根
energy_values = []
wavefunctions = []

for n in range(5):  # 求解前5个状态
    left, right = 0.1 + n, 20 + n
    sol = root_scalar(matching_condition, bracket=[1.5*n, 1.5*n+1.5], method='brentq')
    energy_values.append(sol.root)

    # 使用求得的 E 计算波函数
    sol_wave = solve_ivp(schrodinger, [0, N], [1e-3, 1e-2], args=(sol.root,),
                         t_eval=np.arange(0, N))
    wavefunctions.append(sol_wave.y[0])

# 输出前5个本征能量
print("energy:")
for i, E in enumerate(energy_values):
    print(f"第 {i} 态 能量 E{i}: {E:.6f}")