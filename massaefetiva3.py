import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eig
from scipy.sparse import diags

n = 2000
dx = 1/2000
x = np.linspace(0, 1, n+1)
m = np.linspace(0.1, 1, n+1)

def V(x):
    return 0*x

V = V(x)

# Construindo o vetor A
A = ((1/2)*(1/dx**2)) * (1/m[:-2] + 2/m[1:-1] + 1/m[2:] + V[1:-1])
A = np.insert(A, 0, ((1/2)*(1/dx**2)) * (1/m[0] + 2/m[1] + V[0]))
A = np.append(A, ((1/2)*(1/dx**2)) * (1/m[-2] + 2/m[-1] + V[-1]))

# Construindo os vetores B e C
B = ((-1/2)*(1/dx**2)) * (1/m[1:-1] + 1/m[2:])
B = np.append(B, ((-1/2)*(1/dx**2)) * (1/m[-1]))

C = ((-1/2)*(1/dx**2)) * (1/m[1:-1] + 1/m[:-2])
C = np.insert(C, 0, ((-1/2)*(1/dx**2)) * (1/m[0]))

# Criando a matriz tridiagonal
diagonais = [C, A, B]
offsets = [-1, 0, 1]
matriz = diags(diagonais, offsets).toarray()

# Diagonalizando a matriz
E, psi = eig(matriz)

# Ordenando os autovalores e autofunções
idx = E.argsort()
E = E[idx]
psi = psi[:, idx]

# Normalizando as autofunções
psi = psi / np.linalg.norm(psi, axis=0)

print(f'Autovalores:\n{E}')
print(f'Autovetores:\n{psi}')

# Plotando as autofunções

plt.plot(x, psi[:, 0])
plt.plot(x, psi[:, 1])
plt.plot(x, psi[:, 2])
plt.plot(x, psi[:, 3])
plt.show()
