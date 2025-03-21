from collections.abc import Callable

def solve_3diag(a: float, b: float, c: float, f: list[float]) -> list[float]:
    n = len(f)
    alpha, beta = [0] * n, [0] * n
    alpha[0] = -b / c
    beta[0] = f[0] / c

    for i in range(1, n):
        alpha[i] = -b / (a * alpha[i - 1] + c)
        beta[i] = (f[i] - a * beta[i - 1]) / (a * alpha[i - 1] + c)
    
    x = [0] * n
    x[n - 1] = beta[n - 1]
    for i in range(n - 2, -1, -1):
        x[i] = alpha[i] * x[i + 1] + beta[i]
    return x


def solve_DE_01(func: Callable[[float], float], a: float, b: float, n: int) -> list[float]:
    h = 1 / n
    f = [func(i * h) for i in range(1, n)]
    f[0] += a / h / h
    f[n - 2] += b / h / h

    y = solve_3diag(-1, -1, 2, f)
    for i in range(n - 1):
        y[i] *= h * h
    
    y = [a] + y + [b]
    return y