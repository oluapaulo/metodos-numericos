import numpy as np
import matplotlib.pyplot as plt
def criar_matriz_diagonal(x, dx, pot, kg):
    n = len(x)
    matriz = np.zeros((n, n))

    for i in range(n):
        v1 = pot(x[i])
        m1 = kg(x[i])

        m2 = kg(x[i] + dx)

        m3 = kg(x[i] - dx)

        valor1 = (1/(2*dx**2)) * (1/m2 + 2/m1 + 1/m3 + v1)
        valor2 = -(1/(2*dx**2)) * (1/m1 + 1/m2)
        valor3 = -(1/(2*dx**2)) * (1/m1 + 1/m3)

        matriz[i, i] = valor1

        for i in range(2000):
            if i <= 1999:
                matriz[i, i + 1] = valor2
                matriz[i + 1, i] = valor3

    return matriz


n = 2000
dx = 1/n
x = np.linspace(0.2, 1.2, n+1)

def pot(x):
    return 1000*(x - 7/10)**2

def kg(x):
    return x**3

matriz_resultante = criar_matriz_diagonal(x, dx, pot, kg)

w, m = np.linalg.eig(matriz_resultante)

plt.plot(x, m[:, 0])
#plt.plot(x, m[:, 1])
#plt.plot(x, m[:, 2])
#plt.plot(x, m[:, 3])
plt.show()



