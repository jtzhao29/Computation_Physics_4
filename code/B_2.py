from B_1 import H_matrix,phi,plot_wavefunction
import numpy as np
import matplotlib.pyplot as plt


def calculate_width(phi_wave: np.ndarray, x: np.ndarray) -> float:
    # 计算波包宽度 w(t)
    rho = np.abs(phi_wave)**2  # 粒子密度
    width = np.sqrt(np.sum((x - 100)**2 * rho))  # 计算加权平均距离的平方根
    return width


if __name__ == "__main__":
    time = np.linspace(0, 50, 500)
    width = []
    L = 200
    phi0 = np.zeros(L, dtype=complex)
    phi0[100] = 1.0
    x = np.arange(0,L)
    for t in time:
        phi_t = phi(L,t,phi0)
        width.append(calculate_width(phi_t,x))
    width = np.array(width)
    plt.plot(time,width)    
    plt.xlabel("time",fontsize=20)
    plt.ylabel("width",fontsize=20)
    plt.title("width of wavefunction",fontsize=20)
    plt.grid(True)
    path = f"figure/B_2.png"
    plt.savefig(path)
    plt.show()
