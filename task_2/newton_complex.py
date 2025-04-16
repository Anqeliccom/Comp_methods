import numpy as np
import matplotlib.pyplot as plt

def g(z):
    return z**3 - 1

def g_der(z):
    return 3 * z**2

def newtons_method(f, f_der, x0, delta):
    x1 = x0 - f(x0) / f_der(x0)
    while abs(x1 - x0) > delta:
        x0 = x1
        x1 = x0 - f(x0) / f_der(x0)
    return x1

def classify(z, nodes, delta = 1e-2):
    for ind, node in enumerate(nodes):
        if abs(z - node) < delta:
            return ind
    return -1

N = 800
x = np.linspace(-2, 2, N)
y = np.linspace(-2, 2, N)
X, Y = np.meshgrid(x, y)
Z = X + 1j * Y
image = np.zeros((N, N, 3))

roots = [1, -0.5 + 0.866025j, -0.5 - 0.866025j]
colors = [[1, 0, 0], [0, 1, 0], [0, 0, 1], [0, 0, 0]] # к, з, с, ч

for i in range(N):
    for j in range(N):
        zf = newtons_method(g, g_der, Z[i, j], 1e-3)
        index = classify(zf, roots)
        image[i, j] = colors[index]

plt.figure(figsize=(8, 8))
plt.imshow(image, extent=(-2, 2, -2, 2))
plt.title("Бассейн Ньютона для z^3 - 1")
plt.xlabel("Re")
plt.ylabel("Im")
plt.grid(False)
plt.show()