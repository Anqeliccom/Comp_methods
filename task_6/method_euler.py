import numpy as np
import matplotlib.pyplot as plt

def system(t, y):
    u, v = y
    u_der = 998 * u + 1998 * v
    v_der = -999 * u - 1999 * v
    return np.array([u_der, v_der])

def exact_sol(t):
    return 2 * np.exp(-t) - np.exp(-1000 * t), np.exp(-1000 * t) - np.exp(-t)

def explicit_euler(f, y0, t):
    y = np.zeros((len(t), len(y0)))
    y[0] = y0
    for i in range(len(t) - 1):
        h = t[i + 1] - t[i]
        y[i + 1] = y[i] + h * f(t[i], y[i])
    return y

def implicit_euler(f, y0, t):
    y = np.zeros((len(t), len(y0)))
    y[0] = y0
    for i in range(len(t) - 1):
        h = t[i + 1] - t[i]
        a11 = 1 - 998 * h
        a12 = -1998 * h
        a21 = 999 * h
        a22 = 1 + 1999 * h
        det = a11 * a22 - a12 * a21
        y[i + 1] = [(a22 * y[i, 0] - a12 * y[i, 1]) / det, (-a21 * y[i, 0] + a11 * y[i, 1]) / det]
    return y

def calculate_order(method, f, y0, h):
    t_h = np.arange(0, 0.1, h)
    y_num_h = method(f, y0, t_h)[-1]
    y_exact_h = np.array(exact_sol(t_h[-1]))
    error_h = np.abs(y_num_h - y_exact_h)

    t_h2 = np.arange(0, 0.1, h / 2)
    y_num_h2 = method(f, y0, t_h2)[-1]
    y_exact_h2 = np.array(exact_sol(t_h2[-1]))
    error_h2 = np.abs(y_num_h2 - y_exact_h2)

    print(np.log2(error_h / error_h2))

if __name__ == "__main__":
    y0 = [1.0, 0.0]
    t = np.linspace(0, 0.1, 100)
    u_exact, v_exact = exact_sol(t)

    plt.figure(figsize=(12, 6))

    y_euler = explicit_euler(system, y0, t)
    plt.plot(t, y_euler[:, 0], 'b--', label='Явный Эйлер (u)')
    plt.plot(t, y_euler[:, 1], 'g--', label='Явный Эйлер (v)')

    y_implicit = implicit_euler(system, y0, t)
    plt.plot(t, y_implicit[:, 0], 'r-', label='Неявный Эйлер (u)')
    plt.plot(t, y_implicit[:, 1], 'm-', label='Неявный Эйлер (v)')

    plt.plot(t, u_exact, 'k--', label='Точное (u)')
    plt.plot(t, v_exact, 'k--', label='Точное (v)')

    plt.title('Сравнение методов Эйлера')
    plt.xlabel('Время')
    plt.ylabel('Значения u и v')
    plt.legend()
    plt.grid(True)
    plt.show()

    print("\nЯвный метод Эйлера:")
    calculate_order(explicit_euler, system, y0, 0.001)

    print("\nНеявный метод Эйлера:")
    calculate_order(implicit_euler, system, y0, 0.001)
