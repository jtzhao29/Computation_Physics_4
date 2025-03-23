import numpy as np

def potential(x:float) -> float:
    return 0.5 * x**2 + 0.25 * x**4

def func(phi: np.ndarray, x: float, E: float) -> np.ndarray:
    dphi = phi[1] 
    ddphi = 2 * (potential(x) - E) * phi[0]  
    return np.array([dphi, ddphi])  

def runge_kutta_4(f, phi0: np.ndarray[float], x0: float, xf: float, dx: float, e: float) -> np.ndarray:
    x = x0
    phip = phi0.astype(float)
    length = int((xf - x0) / dx) + 1
    trajectory = np.zeros((length, 3))
    time = 0
    while x < xf:
        phi,dphi = phip[0],phip[1]
        trajectory[time,:] = np.array([x, phi, dphi])
        k1 = dx * f(phip, x, e)
        k2 = dx * f(phip + 0.5 * k1, x + 0.5 * dx, e)
        k3 = dx * f(phip + 0.5 * k2, x + 0.5 * dx, e)
        k4 = dx * f(phip + k3, x + dx, e)
        phip += (k1 + 2 * k2 + 2 * k3 + k4) / 6
        x += dx
        time += 1
    return trajectory


def interative(E0: float, E1: float, L: float, max_iter: int, dx: float, phi0: np.ndarray, tolerate: float):
    x0 = 0
    xf = L
    for i in range(max_iter):
        phi_E0 = runge_kutta_4(func, phi0, x0, xf, dx, E0)[-1, 1]  # 获取 psi(L)
        phi_E1 = runge_kutta_4(func, phi0, x0, xf, dx, E1)[-1, 1]  # 获取 psi(L)
        
        if abs(phi_E1) < tolerate:  
            return E1, runge_kutta_4(func, phi0, x0, xf, dx, E1)
        
        E_new = E1 - phi_E1 * (E1 - E0) / (phi_E1 - phi_E0)
        E0, E1 = E1, E_new

    raise ValueError("割线法不收敛")

def normalize_wavefunction(trajectory):
    x_vals = trajectory[:, 0]
    psi_vals = trajectory[:, 1]
    norm = np.sqrt(np.trapz(psi_vals**2, x_vals))
    return psi_vals / norm  # 返回归一化的波函数

if __name__ == "__main__":
    # 基态，设x=0处波函数为1，导数为0
    phi0 = np.array([1.0,0.0])
    E0 = 0
    E1 = 0.2
    L = 10
    max_iter = 100
    dx = 0.01
    tolerate = 1e-6
    E, trajectory = interative(E0, E1, L, max_iter, dx, phi0, tolerate)
    psi = normalize_wavefunction(trajectory)
    print(f"ground state energy:{E}")
    print(f"param of normal:{np.trapz(psi**2, trajectory[:, 0])}")
