from B_1 import H_matrix,phi,plot_wavefunction
import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':
    time = np.linspace(0, 50, 500)
    rho = []
    L=200
    phi0 = np.zeros(L, dtype=complex)
    phi0[100] = 1.0
    for t in time:
        phi_t = phi(200,t,phi0)
        rho.append(np.abs(phi_t[100])**2)
    rho = np.array(rho)
    plt.plot(time,rho)
    plt.xlabel("time",fontsize=20)
    plt.ylabel("density of center",fontsize=20)
    plt.title("density of wavefunction in center",fontsize=20)
    plt.grid(True)
    path = f"figure/B_3.png"
    plt.savefig(path)
    plt.show()
