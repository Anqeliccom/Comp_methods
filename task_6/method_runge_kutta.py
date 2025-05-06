import numpy as np
import matplotlib.pyplot as plt

def test_equation(t, y):
    return -3 * y * t ** 2

def exact_sol(t):
    return np.exp(-t ** 3)

def lorenz_attractor(t, vector, sigma, b, r):
    x, y, z = vector
    x_der = sigma * (y - x)
    y_der = x * (r - z) - y
    z_der = x * y - b * z
    return np.array([x_der, y_der, z_der])

def runge_kutta_4(f, y0, t, args=()):
    y = np.zeros((len(t), len(y0)))
    y[0] = y0

    for i in range(len(t) - 1):
        h = t[i + 1] - t[i]
        k1 = f(t[i], y[i], *args)
        k2 = f(t[i] + h / 2, y[i] + h * k1 / 2, *args)
        k3 = f(t[i] + h / 2, y[i] + h * k2 / 2, *args)
        k4 = f(t[i] + h, y[i] + h * k3, *args)
        y[i + 1] = y[i] + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
    return y

def calculate_order(h):
    y0 = [1.0]

    t_h = np.arange(0, 2, h)
    y_num_h = runge_kutta_4(test_equation, y0, t_h)[:, 0]
    y_exact_h = exact_sol(t_h)
    error_h = np.max(np.abs(y_num_h - y_exact_h))

    t_h2 = np.arange(0, 2, h / 2)
    y_num_h2 = runge_kutta_4(test_equation, y0, t_h2)[:, 0]
    y_exact_h2 = exact_sol(t_h2)
    error_h2 = np.max(np.abs(y_num_h2 - y_exact_h2))

    print(np.log2(error_h / error_h2))

if __name__ == "__main__":
    sigma = 10.0
    b = 8.0 / 3.0
    r = 25.0
    initial_cond = np.array([1.0, 1.0, 1.0])
    t = np.arange(0, 50, 0.01)
    solution = runge_kutta_4(lorenz_attractor, initial_cond, t, args=(sigma, b, r))

    fig = plt.figure(figsize=(12, 9))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(solution[:, 0], solution[:, 1], solution[:, 2], lw=0.5)
    ax.set_title("Аттрактор Лоренца")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    plt.show()
    calculate_order(0.01)
