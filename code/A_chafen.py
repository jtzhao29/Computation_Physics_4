import numpy as np
import matplotlib.pyplot as plt

def matrix(L:float,dx:float)->tuple[np.ndarray,np.ndarray]:
    N = int(2*L/dx)
    x = np.linspace(-L,L,N)
    # 势能项加到对角线上
    V = 0.5 * x**2 + 0.1 * x**4
    H = np.zeros((N, N))
    for i in range(N):
        if i > 0:
            H[i, i-1] = -1 / (2 * dx**2)
        H[i, i] = V[i] + 1 / (dx**2)
        if i < N - 1:
            H[i, i+1] = -1 / (2 * dx**2)
    return H,x
def matrix_jianxie(L:float,dx:float)->tuple[np.ndarray,np.ndarray]:
    N = int(2*L/dx)
    x = np.linspace(-L,L,N)
    # 势能项加到对角线上
    V = 0.5 * x**2 
    H = np.zeros((N, N))
    for i in range(N):
        if i > 0:
            H[i, i-1] = -1 / (2 * dx**2)
        H[i, i] = V[i] + 1 / (dx**2)
        if i < N - 1:
            H[i, i+1] = -1 / (2 * dx**2)
    return H,x

def eigen_and_norm(H:np.ndarray, x:np.ndarray)->tuple[np.ndarray,np.ndarray]:
    eigvals, eigvecs = np.linalg.eigh(H)    
    for n in range(eigvecs.shape[1]):
        norm = np.sqrt(np.trapz(eigvecs[:, n]**2, x))
        eigvecs[:, n] /= norm
    eigvals = np.round(eigvals, 4)  # 精确到1e-4
    return eigvals, eigvecs

def plot_wavefunction_abs(x:np.ndarray,wavefunction:np.ndarray,i:int)->None:
    plt.figure(figsize=(10, 6))
    plt.plot(x,np.abs(wavefunction), label=f"|ψ(x)|")
    plt.xlabel("x",fontsize=20)
    plt.ylabel("Amplitude",fontsize=20)
    plt.title(f"wavefunction |ψ(x)| (Normalized) of energy level {i}",fontsize=20)
    plt.legend(loc='upper right')
    plt.grid()
    
    path = f"figure/A_chafen_{i}.png"
    plt.savefig(path)  
    plt.show()  

def plot_wavefunction_compare(x:np.ndarray,wavefunction:np.ndarray,wavefunction1:np.ndarray,i:int)->None:
    plt.figure(figsize=(10, 6))
    plt.plot(x,wavefunction, label=f"|ψ(x)| of question A")
    plt.plot(x,wavefunction1, label=f"|ψ(x)| of harmonic oscillator")
    plt.xlabel("x",fontsize=20)
    plt.ylabel("Amplitude",fontsize=20)
    plt.title(f"wavefunction ψ(x) (Normalized) of energy level {i}",fontsize=20)
    plt.legend(loc='upper right')
    plt.grid()
    
    path = f"figure/A_chafen_compare_{i}.png"
    plt.savefig(path)  
    plt.show()  


if __name__ == "__main__":
    L = 10
    dx = 0.01
    H ,x= matrix(L,dx)
    eigvals,eigvecs = eigen_and_norm(H,x)
    # ground_state_energy = eigvals[0]
    # ground_state_wavefunction = eigvecs[:, 0]
    # plot_wavefunction_abs(x,ground_state_wavefunction,0)
    # print(f"Ground State Energy: {ground_state_energy}")
    # for i in range(1, 5):
    #     energy = eigvals[i]
    #     wavefunction = eigvecs[:, i]
    #     print(f"Energy Level {i}: {energy}")
    #     plot_wavefunction_abs(x,wavefunction,i)
    H1,x = matrix_jianxie(L,dx)
    eigvals1,eigvecs1 = eigen_and_norm(H1,x)
    for i in range(0, 5):
        energy = eigvals[i]
        energy1 = eigvals1[i]
        wavefunction0 = eigvecs[:, i]
        wavefunction1 = eigvecs1[:, i]
        print(f"Energy Level {i} A:{energy}")
        print(f"Energy Level {i} of harmonic oscillator:{energy1}")
        plot_wavefunction_compare(x,wavefunction0,wavefunction1,i)