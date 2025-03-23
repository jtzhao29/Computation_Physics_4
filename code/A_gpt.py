import numpy as np
import matplotlib.pyplot as plt

# 定义势能函数 V(x)
def potential(x: float) -> float:
    return 0.5 * x**2 + 0.25 * x**4

# 定义 φ 的微分方程组（薛定谔方程）
def func(phi: np.ndarray, x: float, E: float) -> np.ndarray:
    dphi = phi[1]  # ψ'(x)
    ddphi = 2 * (potential(x) - E) * phi[0]  # ψ''(x)
    return np.array([dphi, ddphi])  # 返回 [ψ'(x), ψ''(x)]

# 实现四阶龙格-库塔法求解
def runge_kutta_4(f, phi0: np.ndarray, x0: float, xf: float, dx: float, e: float) -> np.ndarray:
    x = x0
    phip = phi0.astype(float)  # 波函数初始条件
    length = int((xf - x0) / dx) + 1  # 计算步数
    trajectory = np.zeros((length, 3))  # 记录解轨迹
    
    for time in range(length):  # 使用步进循环而非 while，避免越界
        trajectory[time, :] = np.array([x, phip[0], phip[1]])  # 记录 x, ψ(x), ψ'(x)
        k1 = dx * f(phip, x, e)
        k2 = dx * f(phip + 0.5 * k1, x + 0.5 * dx, e)
        k3 = dx * f(phip + 0.5 * k2, x + 0.5 * dx, e)
        k4 = dx * f(phip + k3, x + dx, e)
        phip += (k1 + 2 * k2 + 2 * k3 + k4) / 6
        x += dx  # 更新当前 x 值
    return trajectory  # 返回轨迹

# 割线法，寻找满足边界的本征能量
def interative(E0: float, E1: float, L: float, max_iter: int, dx: float, phi0: np.ndarray, tolerate: float):
    x0 = 0
    xf = L
    for i in range(max_iter):
        # 计算 f(E0) 和 f(E1)，即边界点 ψ(L) 的值
        phi_E0 = runge_kutta_4(func, phi0, x0, xf, dx, E0)[-1, 1]  # ψ(L) for E0
        phi_E1 = runge_kutta_4(func, phi0, x0, xf, dx, E1)[-1, 1]  # ψ(L) for E1
        
        # 检查是否已收敛
        if abs(phi_E1) < tolerate:
            return E1, runge_kutta_4(func, phi0, x0, xf, dx, E1)  # 返回能量及波函数轨迹
        
        # 割线法更新公式
        E_new = E1 - phi_E1 * (E1 - E0) / (phi_E1 - phi_E0)
        E0, E1 = E1, E_new  # 更新能量值
    
    raise ValueError(f"割线法不收敛。最后的 E0={E0}, E1={E1}, 误差={phi_E1}")  # 调试信息

# 对波函数进行归一化
def normalize_wavefunction(trajectory):
    x_vals = trajectory[:, 0]  # x 的值
    psi_vals = trajectory[:, 1]  # ψ(x) 的值
    norm = np.sqrt(np.trapz(psi_vals**2, x_vals))  # 通过梯形积分计算归一化因子
    return psi_vals / norm  # 返回归一化的波函数

if __name__ == "__main__":
    # 设置基态解的初始条件和参数
    phi0 = np.array([1.0, 0.0])  # 偶函数波函数初始条件：ψ(0)=1, ψ'(0)=0
    E0 = 0  # 初始猜测能量
    E1 = 0.2  # 第二个猜测能量
    L = 10  # 数值积分的边界（x 的范围）
    max_iter = 100  # 最大迭代次数
    dx = 0.01  # 数值积分步长
    tolerate = 1e-6  # 收敛条件
    
    # 调用割线法求解基态能量
    try:
        E, trajectory = interative(E0, E1, L, max_iter, dx, phi0, tolerate)
        psi = normalize_wavefunction(trajectory)  # 对波函数归一化
        print(f"Ground state energy: {E:.6f}")  # 打印基态能量
        print(f"Normalization check (should be ~1): {np.trapz(psi**2, trajectory[:, 0]):.6f}")
    except ValueError as e:
        print(e)

    # 绘制归一化后的波函数 |ψ(x)|
    x_vals = trajectory[:, 0]
    plt.plot(x_vals, psi, label="Ground State ψ(x)")
    plt.xlabel("x")
    plt.ylabel("ψ(x)")
    plt.title(f"Ground state wavefunction, E={E:.6f}")
    plt.legend()
    plt.grid()
    plt.show()