from DE_solver import solve_DE_01
from collections.abc import Callable
import math


def func(x: float) -> float:
    return math.sin(x)


def right_f(x: float) -> float:
    return math.sin(x)


def C_norm_01(func: Callable[[float], float], discrete_func: list[float]):
    n = len(discrete_func) - 1
    h = 1 / (n)
    return max([abs(func(h * i) - discrete_func[i]) for i in range(n + 1)])


def L2_norm_01(func: Callable[[float], float], discrete_func: list[float]):
    n = len(discrete_func) - 1
    h = 1 / (n)
    norm = 0
    for i in range(n + 1):
        norm += (func(h * i) - discrete_func[i]) ** 2
    return math.sqrt(norm)


def test(n: int):
    res = solve_DE_01(right_f, func(0), func(1), n)
    return C_norm_01(func, res), L2_norm_01(func, res)


print(f" {'N':<6} | {'C-norm':>25} | {'L2-norm':>25}")
for n in range(3, 20):
    print(" {:<6} | {:>25} | {:>25}".format(2 ** n, *test(2 ** n)))