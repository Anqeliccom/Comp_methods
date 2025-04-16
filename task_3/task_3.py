import math
import matplotlib.pyplot as plt

def f(x): return 1 / (1 + 25 * x ** 2)

def xi(n, a, b):
    return [2 * i / n - 1 for i in range(n + 1)]

def chebyshev_xi(n, a, b):
    return [0.5 * (a + b) + 0.5 * (b - a) * math.cos(math.pi * (2 * k + 1) / (2 * n + 2)) for k in range(n + 1)]

def div_diff(x, y):
    n = len(x)
    table = [[0.0] * n for _ in range(n)]
    for i in range(n):
        table[i][0] = y[i]

    for j in range(1, n):
        for i in range(n - j):
            table[i][j] = (table[i + 1][j - 1] - table[i][j - 1]) / (x[i + j] - x[i])
    return [table[0][j] for j in range(n)]

def newtons_polynomial(x, y):
    d = div_diff(x, y)
    def P(xp):
        res, product = d[0], 1
        for i in range(1, len(x)):
            product *= (xp - x[i - 1])
            res += d[i] * product
        return res
    return P

def graph(ns, node_func, title, a, b, h):
    xp = [a + i * h for i in range(int((b - a) / h) + 1)]
    max_errors = []

    for n in ns:
        xn = node_func(n, a, b)
        yn = [f(x) for x in xn]
        yp = [abs(newtons_polynomial(xn, yn)(x) - f(x)) for x in xp]
        max_errors.append(max(yp))

    plt.figure(figsize=(10, 6))
    plt.plot(ns, max_errors, marker='o')
    plt.title(title)
    plt.xlabel("n")
    plt.ylabel("Максимальная погрешность")
    plt.grid()
    plt.show()

nl = list(range(3, 16))
a, b = -1, 1
h = 0.001
graph(nl, xi, "Заданные узлы", a, b, h)
graph(nl, chebyshev_xi, "Узлы Чебышёва", a, b, h)