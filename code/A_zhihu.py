import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import spdiags
from scipy.sparse.linalg import eigsh

# Step 1: 设置基本参数
X = 10          # 区间长度
N = 200         # 格点数
dx = 2 * X / N  # 空间步长
x = np.linspace(-X, X, N)

n = 2           # 势函数次幂
V = 0.5 * x**n +0.25*x**4 # 势场

# Step 2: 构造哈密顿算符
diagonal_V = spdiags(V, 0, N, N)

# 构造哈密顿算符的三对角形式
H = np.zeros((N, 3))
H[:, 0] = -0.5 / dx**2    # 次对角线（下）
H[:, 1] = 1 / dx**2       # 主对角线
H[:, 2] = -0.5 / dx**2    # 次对角线（上）

# 生成稀疏矩阵
H_sparse = spdiags(H, [-1, 0, 1], N, N)

# 最终的哈密顿算符矩阵
Hamiltonian = diagonal_V + H_sparse

# Step 3: 求解本征值和本征函数
E = 9  # 要求解的能级数
eigenvalues, eigenvectors = eigsh(Hamiltonian, k=E, sigma=0)  # sigma=0 表示寻找最小的本征值

# Step 4: 绘图
plt.figure(figsize=(10, 6))

for i in range(E):
    plt.plot(x, eigenvectors[:, i]**2, label=f'Energy Level {i} = {eigenvalues[i]:.4f}')

plt.title("Wave Functions for Bound States")
plt.xlabel("x")
plt.ylabel("ψ(x)")
plt.legend()
plt.grid(True)
plt.show()

print(f"Eigenvalues: {eigenvalues}")

