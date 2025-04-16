import math

def f1(x):
    if x == 0:
        return math.pi
    if x == 1:
        return 5 * math.pi
    return math.sin(math.pi * x ** 5) / (x ** 5 * (1 - x))

def f2(x):
    if x == float('inf'):
        return 0
    return math.exp(-math.sqrt(x) + math.sin(x / 10))

def make_replace(f):
    def new_f(t):
        if t == 1:
            return f(float('inf'))
        return f(t / (1 - t)) / ((1 - t) ** 2)
    return new_f

def simpson_method(f, a, b, n):
    if n % 2 != 0:
        n += 1
    if b == float("inf"):
        f = make_replace(f)
        b = 1
    h = (b - a) / n
    x = [a + i * h for i in range(n + 1)]
    fx = [f(xi) for xi in x]

    sum_odd  = sum(fx[i] for i in range(1, n + 1, 2))
    sum_even = sum(fx[i] for i in range(2, n, 2))

    result = fx[0] + 2 * sum_even + 4 * sum_odd + fx[n]
    return result * h / 3

def simpson_with_accurate(f, a, b, delta, exact_value):
    n = 2
    curr = simpson_method(f, a, b, n)
    while abs(curr - exact_value) >= delta:
        n *= 2
        curr = simpson_method(f, a, b, n)
    return curr, n

def calculate_order(f, a, b, n):
    integral_n    = simpson_method(f, a, b, n)         # I(h)
    integral_2n   = simpson_method(f, a, b, 2 * n)     # I(h/2)
    integral_4n   = simpson_method(f, a, b, 4 * n)     # I(h/4)

    p = math.log2(abs(integral_2n - integral_n)  / abs(integral_4n - integral_2n) )
    return p

a1, b1 = 0, 1
a2, b2 = 0, float("inf")
e = 1e-8
exact_value1 = 8.03491057
exact_value2 = 2.98100345
result1, n1 = simpson_with_accurate(f1, a1, b1, e, exact_value1)
print(f"Интеграл f1 на [{a1}, {b1}] = {result1:.10f}, при n = {n1}")

result2, n2 = simpson_with_accurate(f2, a2, b2, e, exact_value2)
print(f"Интеграл f2 на [{a2}, {b2}] = {result2:.10f}, при n = {n2}")

p1 = calculate_order(f1, a1, b1, n1)
print(f"Порядок аппроксимации: p = {p1:.2f}")

p2 = calculate_order(f2, a2, b2, n2)
print(f"Порядок аппроксимации: p = {p2:.2f}")