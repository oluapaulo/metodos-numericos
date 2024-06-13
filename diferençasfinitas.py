import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh_tridiagonal

n = 2000
dx = 1/n
x = np.linspace(0, 1, n+1)

def pot(x):
    return 1000*(x - 1/2)**2

v = pot(x)



d = 2/dx**2 + v[1:-1]
e = -1/dx**2 * np.ones(len(d) - 1)
w, m = eigh_tridiagonal(d, e)
plt.plot(x[1:-1], m[:, 0])
plt.plot(x[1:-1], m[:, 1])
plt.plot(x[1:-1], m[:, 2])
plt.plot(x[1:-1], m[:, 3])
plt.show()

plt.bar(np.arange(0,10,1), w[0:10])
#plt.show()