import math

e = 0.001

def f(x):
    return math.tan(x) - x

def f_der(x):
    return (1 / math.cos(x))**2 - 1

def g(z):
    return z**3 - 1

def g_der(z):
    return 3 * z**2

def bisection_method(f, a, b, delta, count = 0):
    count += 1
    c = (a + b) / 2
    if abs(b - a) <= 2 * delta or f(c) == 0:
        return c
    if f(a) * f(c) < 0:
        return bisection_method(f, a, c, delta, count)
    return bisection_method(f, c, b, delta, count)

def simple_iteration_method(a, b, delta, n):
    x0 = (a + b) / 2
    x1 = math.atan(x0)
    while abs(x1 - x0) > delta:
        x0 = x1
        x1 = math.atan(x0) + math.pi * n
    return x1

def newtons_method(f, f_der, x0, delta):
    x1 =  x0 - f(x0) / f_der(x0)
    while abs(x1 - x0) > delta:
        x0 = x1
        x1 = x0 - f(x0) / f_der(x0)
    return x1

def secant_method(f, x0, x1, delta):
    for _ in range(3): # чтобы развести точки
        x2 = x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0))
        x0, x1 = x1, x2

    x2 = x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0))
    while abs(x1 - x0) > delta:
        x0, x1 = x1, x2
        x2 = x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0))
    return x2

def tests():
    print("Bisection method")
    print("-" * 20)
    for i in range(1, 11):
        a = (i - 0.5) * math.pi + e
        b = (i + 0.5) * math.pi - e
        print(bisection_method(f, a, b, 0.0001))

    print("\nSimple iteration method")
    print("-" * 20)
    for i in range(1, 11):
        a = (i - 0.5) * math.pi + e
        b = (i + 0.5) * math.pi - e
        print(simple_iteration_method(a, b, 0.0001, i))

    print("\nNewton's method")
    print("-" * 20)
    for i in range(1, 11):
        x0 = (i + 0.5) * math.pi - e
        print(newtons_method(f, f_der, x0, 0.00001))

    print("\nSecant method")
    print("-" * 20)
    for i in range(1, 11):
        a = (i + 0.5) * math.pi - 2 * e
        b = (i + 0.5) * math.pi - e
        print(secant_method(f, a, b, 0.0001))

    print("\nNewton's method (z^3-1)")
    print("-" * 20)
    print(newtons_method(g, g_der, 1.5, 0.001))
    print(newtons_method(g, g_der, 0.5 + 0.8j, 0.001))
    print(newtons_method(g, g_der, 0.5 - 0.8j, 0.001))

tests()