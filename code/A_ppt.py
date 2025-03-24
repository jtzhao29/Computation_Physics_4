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



def normalize_wavefunction(trajectory):
    x_vals = trajectory[:, 0]
    psi_vals = trajectory[:, 1]
    norm = np.sqrt(np.trapz(psi_vals**2, x_vals))
    return psi_vals / norm  # 返回归一化的波函数

def fff(E:float,L:float,dx:float,phi0:np.ndarray)->float:
    X01 = -L
    X02 = L
    xf =0.0
    trajectory1 = runge_kutta_4(func, phi0, X01, xf, dx, E)
    trajectory2 = runge_kutta_4(func, phi0, X02, xf, -dx, E)
    return (trajectory1[-1, 1] - trajectory2[-1, 1])**2 + (trajectory1[-1, 2] - trajectory2[-1, 2])**2

def interative(E0: float, E1: float, L: float, max_iter: int, dx: float, phi0: np.ndarray, tolerate: float)->float:
    E_new = E1-fff(E1,L,dx,phi0)*(E1-E0)/(fff(E1,L,dx,phi0)-fff(E0,L,dx,phi0))
    for i in range(max_iter):
        if abs(fff(E_new,L,dx,phi0)) < tolerate:
            return E_new
        E0,E1 = E1,E_new
        E_new = E1-fff(E1,L,dx,phi0)*(E1-E0)/(fff(E1,L,dx,phi0)-fff(E0,L,dx,phi0))
    raise ValueError("111")

if __name__ == "__main__":
    # 设+-L处波函数为0，导数为0
    phi0 = np.array([0.0,0.0])
    E0 = 0.0
    E1 = 1.0
    L = 90
    max_iter = 100
    dx = 0.01
    tolerate = 1e-6
    E = interative(E0, E1, L, max_iter, dx, phi0, tolerate)
    print(E)
    
