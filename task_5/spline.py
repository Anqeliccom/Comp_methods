import numpy as np
import matplotlib.pyplot as plt
import task_5

def func(x):
    return 1 / (1 + 25 * x ** 2)

def get_coefficients(x, y, lcond, rcond):
    n = len(x) - 1
    h = [x[i + 1] - x[i] for i in range(n)]

    A = []
    B = []
    C = []
    rhs = []

    for i in range(1, n):
        A.append(h[i - 1])
        C.append(2 * (h[i - 1] + h[i]))
        B.append(h[i])
        rhs.append(3 * ((y[i + 1] - y[i]) / h[i] - (y[i] - y[i - 1]) / h[i - 1]))

    if lcond[0]:  # первая производная
        C = [1.0] + C
        B = [-0.5] + B
        rhs = [(1.5 / h[0]) * ((y[1] - y[0]) / h[0] - lcond[1])] + rhs
    else:  # вторая производная
        C = [1.0] + C
        B = [0.0] + B
        rhs = [lcond[1]] + rhs

    if rcond[0]:
        A = A + [-2.0]
        C = C + [1.0]
        rhs = rhs + [(3 / h[n-1]) * ((y[n] - y[n - 1]) / h[n - 1] - rcond[1])]
    else:
        A = A + [0.0]
        C = C + [1.0]
        rhs = rhs + [rcond[1]]

    c = task_5.tridiagonal_method(A, C, B, rhs)

    b = [(y[i + 1] - y[i]) / h[i] - h[i] * (2 * c[i] + c[i + 1]) / 3 for i in range(n)]
    d = [(c[i + 1] - c[i]) / (3 * h[i]) for i in range(n)]

    return c, b, d

def get_spline(x, y, c, b, d, xp):
    for i in range(len(x) - 1):
        if x[i] <= xp <= x[i + 1]:
            dx = xp - x[i]
            return y[i] + b[i] * dx + c[i] * dx ** 2 + d[i] * dx ** 3
    return None

def plot_spline_vs_function(f, a, b, n, lcond, rcond):
    x = np.linspace(a, b, n + 1)
    y = f(x)
    c, b_coef, d = get_coefficients(x, y, lcond, rcond)
    xp = np.linspace(a, b, 500)
    sp = [get_spline(x, y, c, b_coef, d, xi) for xi in xp]

    plt.figure(figsize=(10, 6))
    plt.plot(xp, f(xp), label='f(x)')
    plt.plot(xp, sp, label='Spline', linestyle='--')
    plt.legend()
    plt.title('Сравнение функции и сплайна')
    plt.grid()
    plt.show()

def plot_max_errors(f, a, b, lcond, rcond, max_nodes):
    ns = range(10, max_nodes + 1, 10)
    max_errors = []

    for n in ns:
        x = np.linspace(a, b, n + 1)
        y = f(x)
        c, b_coef, d = get_coefficients(x, y, lcond, rcond)
        xp = np.linspace(a, b, 500)
        sp = [get_spline(x, y, c, b_coef, d, xi) for xi in xp]
        max_errors.append(max(abs(si - fxi) for si, fxi in zip(sp, f(xp))))

    plt.figure(figsize=(10, 6))
    plt.plot(ns, max_errors, 'o-')
    plt.xlabel('Количество узлов')
    plt.ylabel('Максимальная погрешность')
    plt.title('Зависимость погрешности от количества узлов')
    plt.grid()
    plt.show()

a, b = -1, 1
lcond = (False, 0.0)  # S''(a) = 0
rcond = (False, 0.0)  # S''(b) = 0

plot_spline_vs_function(func, a, b, 50, lcond, rcond)
plot_max_errors(func, a, b, lcond, rcond, 100)