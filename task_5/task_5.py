import matplotlib.pyplot as plt
import numpy as np

def equation_rhs(x):
    return np.cos(x)

def tridiagonal_method(A, C, B, F):
    k1, m1 = B[0], F[0]
    k2, m2 = A[-1], F[-1]
    alpha, beta = [k1], [m1]
    n = len(F)
    for i in range(n - 1):
        alpha.append(B[i] / (C[i] - A[i] * alpha[i]))
        beta.append((A[i] * beta[i] + F[i]) / (C[i] - A[i] * alpha[i]))

    y = [0.0] * n
    y[-1] = (m2 + k2 * beta[n - 1]) / (1 - k2 * alpha[n - 1])
    for i in range(n - 1, 0, -1):
        y[i - 1] = alpha[i] * y[i] + beta[i]
    return y

def get_matrices(func, lcond, rcond, a, b, n):
    h = (b - a) / n
    X = [a + i * h for i in range(1,n)]

    A = [1 / h ** 2] * (n)
    C = [2 / h ** 2] * (n)
    B = [1 / h ** 2] * (n)
    F = [-func(x) for x in X]
    X = [a] + X + [b]

    if lcond[0]:  # производная
        C = [1.0] + C
        B = [1.0] + B
        F = [-h * lcond[1]] + F
    else:  # функция
        C = [1.0] + C
        B = [0.0] + B
        F = [lcond[1]] + F

    if rcond[0]:
        A = A + [1.0]
        C = C + [1.0]
        F = F + [h * rcond[1]]
    else:
        A = A + [0.0]
        C = C + [1.0]
        F = F + [rcond[1]]

    return A, C, B, F, X

def examples(example_num):
    ex = {
        1: ((False, 0.0), (False, 1.0), lambda x: -np.cos(x) + x / np.pi + 1 / 2), # y(a)=0, y(b)=1
        2: ((True, 0.0), (False, 1.0), lambda x: -np.cos(x) + x - np.pi / 2 + 1), # y'(a)=0, y(b)=1
        3: ((False, 0.0), (True, 0.0), lambda x: -np.cos(x) - x - np.pi / 2), # y(a)=0, y'(b)=0
    }
    return ex.get(example_num, ex[1])

def plot_max_errors(a, b, lcond, rcond, exact_sol, max_nodes):
    ns = range(10, max_nodes + 1, 10)
    max_errors = []

    for n in ns:
        A, C, B, F, X = get_matrices(equation_rhs, lcond, rcond, a, b, n)
        y_num = tridiagonal_method(A, C, B, F)
        y_exact = [exact_sol(x) for x in X]
        max_errors.append(max(abs(num - exact) for num, exact in zip(y_num, y_exact)))

    plt.figure(figsize=(10, 6))
    plt.plot(ns, max_errors, 'o-')
    plt.xlabel('Количество узлов')
    plt.ylabel('Максимальная погрешность')
    plt.title('Зависимость погрешности от количества узлов')
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    a, b = -np.pi / 2, np.pi / 2
    lcond, rcond, exact_solution = examples(1)

    n = 100
    A, C, B, F, X = get_matrices(equation_rhs, lcond, rcond, a, b, n)
    y_num = tridiagonal_method(A, C, B, F)
    y_exact = [exact_solution(x) for x in X]

    plt.figure(figsize=(10, 6))
    plt.plot(X, y_num, label="Численное решение", marker='o', markersize=3, linestyle='-')
    plt.plot(X, y_exact, label="Точное решение")
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Сравнение численного и точного решений')
    plt.legend()
    plt.grid()
    plt.show()

    plot_max_errors(a, b, lcond, rcond, exact_solution, 100)
