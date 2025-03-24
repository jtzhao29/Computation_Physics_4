import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags
from scipy.sparse.linalg import expm
from scipy.sparse import csr_matrix
def H_matrix(L: int) -> np.ndarray:
    # 构造哈密顿量
    H = np.zeros((L, L))
    for i in range(L):
        if i > 0:
            H[i, i-1] = -1
        if i < L - 1:
            H[i, i+1] = -1
    # 周期边界条件
    H[0, L-1] = H[L-1, 0] = -1
    return H

def phi(L: int, t: float, phi0: np.ndarray) -> np.ndarray:
    H = csr_matrix(H_matrix(L)) 
    U_t = expm(-1j * t * H)    
    return U_t @ phi0

def plot_wavefunction(L:int,t:float,phi_wave:np.ndarray)->None:
    # 绘制波函数
    x = np.arange(0,L)
    plt.plot(x, np.abs(phi_wave)**2)
    plt.xlabel(r'$x$', fontsize=20) 
    plt.ylabel(r'$|\varphi(x)|^2$', fontsize=20)  # 注意 `\varphi` 更接近物理学中常用的波函数符号
    plt.title(r'$|\varphi(x)|^2$ at $t = %.1f$' % t, fontsize=20)
    plt.grid(True)
    path = f"figure/B_1_{t}.png"
    plt.savefig(path)
    plt.show()

if __name__ == "__main__":
    time = [90.0]
    L=200
    phi0 = np.zeros(L, dtype=complex)
    phi0[100] = 1.0
    for t in time:
        phi_t = phi(L,t,phi0)
        plot_wavefunction(L,t,phi_t)

