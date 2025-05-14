import math
import cmath
import matplotlib.pyplot as plt

def f(z):
    return z ** 3 - 1

def f_der(z):
    return 3 * z ** 2

def trapezoid_method(f, z):
    z_start = z[:-1]
    z_end = z[1:]
    integral = 0
    for zs, ze in zip(z_start, z_end):
        integral += 0.5 * (f(zs) + f(ze)) * (ze - zs)
    return integral

def integrate_rectangle(f, f_der, p, x0, y0, x1, y1, n):
    def formula(z):
        return (z ** p) * f_der(z) / f(z)

    edges = []
    xs = [x0 + i * (x1 - x0) / n for i in range(n + 1)]
    ys = [y0 + i * (y1 - y0) / n for i in range(n + 1)]

    # x0 -> x1
    edges.append([x + 1j * y0 for x in xs])
    # y0 -> y1
    edges.append([x1 + 1j * y for y in ys])
    # x1 -> x0
    edges.append([x + 1j * y1 for x in reversed(xs)])
    # y1 -> y0
    edges.append([x0 + 1j * y for y in reversed(ys)])

    total_integral = 0
    for edge in edges:
        total_integral += trapezoid_method(formula, edge)
    return total_integral / (2j * math.pi)

def calculate_sigma(f, f_der, x0, y0, x1, y1, n):
    N = round(integrate_rectangle(f, f_der, 0, x0, y0, x1, y1, n).real)
    s = [integrate_rectangle(f, f_der, i + 1, x0, y0, x1, y1, n) for i in range(N)]
    sigma = [0] * N
    sigma[0] = -s[0]

    for i in range(1, N):
        sum_ = sum(sigma[i - j - 1] * s[j] for j in range(i))
        sigma[i] = -(sum_ + s[i]) / (i + 1)
    return [c.real for c in sigma]

def create_companion_matrix(c):
    dim = len(c)
    C = [[0 for _ in range(dim)] for _ in range(dim)]
    for i in range(1, dim):
        C[i][i - 1] = 1
    for i in range(dim):
        C[i][dim - 1] = -c[dim - 1 - i]
    return C

def hessenberg_qr(H):
    n = len(H)
    Q = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    R = [row.copy() for row in H]

    for k in range(n - 1):
        x, y = R[k][k], R[k + 1][k]
        r = math.sqrt(x**2 + y**2)
        c = x / r
        s = -y / r

        for j in range(n):
            a, b = Q[j][k], Q[j][k + 1]
            Q[j][k] = c * a - s * b
            Q[j][k + 1] = s * a + c * b

        for j in range(k, n):
            a, b = R[k][j], R[k + 1][j]
            R[k][j] = c * a - s * b
            R[k + 1][j] = s * a + c * b

    return Q, R

def wilkinson_shift(a, b, c, d):
    delta = (a - d) / 2
    if delta == 0:
        return d - math.sqrt(abs(b * c))  # am - |bm|
    denom = delta + math.copysign(abs(cmath.sqrt(delta ** 2 + b * c)), delta)
    return d - (b * c) / denom

def qr_find_eigenvalues(H):
    A = [row.copy() for row in H]
    size = len(A)
    result = []

    idx = size - 1
    while idx >= 0:
        if idx == 0 or abs(A[idx][idx - 1]) < 1e-12 * (abs(A[idx - 1][idx - 1]) + abs(A[idx][idx])): # |bm| < ε(|am−1| + |am|)
            result.append(A[idx][idx])
            idx -= 1
            continue

        a, b = A[idx - 1][idx - 1], A[idx - 1][idx]
        c, d = A[idx][idx - 1], A[idx][idx]

        # λ² - tr*λ + det = 0
        if idx == 1 or abs(A[idx - 1][idx - 2]) < 1e-12:
            tr = a + d
            det = a * d - b * c
            discriminant = tr ** 2 - 4 * det
            root = cmath.sqrt(discriminant)
            result.extend([(tr + root) / 2, (tr - root) / 2])
            idx -= 2
            continue

        shift = wilkinson_shift(a, b, c, d)
        shifted = [[A[i][j] - (shift if i == j else 0) for j in range(idx + 1)] for i in range(idx + 1)] # H_shifted = H - s*I
        Q, R = hessenberg_qr(shifted)

        for i in range(idx + 1):
            for j in range(idx + 1):
                A[i][j] = sum(R[i][k] * Q[k][j] for k in range(idx + 1)) + (shift if i == j else 0) # H = R*Q + s*I

    return result

def find_roots_and_plot(f, f_der, xmin, ymin, xmax, ymax, points):
    poly_coeffs = calculate_sigma(f, f_der, xmin, ymin, xmax, ymax, points)
    comp_matrix = create_companion_matrix(poly_coeffs)
    roots = qr_find_eigenvalues(comp_matrix)

    print("Найденные корни:")
    for root in roots:
        print(root)

    plt.figure(figsize=(6, 6))
    plt.scatter([root.real for root in roots], [root.imag for root in roots], c='red', marker='x')
    plt.grid(True)
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.title("Корни функции в комплексной плоскости")
    plt.show()

if __name__ == "__main__":
    find_roots_and_plot(f, f_der, -10, -10, 10, 10, 500000)