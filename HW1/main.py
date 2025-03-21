from DE_solver import solve_DE_01
from collections.abc import Callable
import matplotlib.pyplot as plt
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
n_list = []
c_norm = []
l2_norm = []
for pow in range(3, 20):
    c, l2 = test(2 ** pow)
    print(" {:<6} | {:>25} | {:>25}".format(2 ** pow, c, l2))
    n_list.append(2 ** pow)
    c_norm.append(c)
    l2_norm.append(l2)

plt.plot(n_list, c_norm, label="С-норма")
plt.plot(n_list, l2_norm, label="L2-норма")
plt.title("C-норма и L2-норма ошибки в зависимости от N")
plt.xlabel("N")
plt.ylabel("норма ошибки")
plt.xscale("log")
plt.yscale("log")
plt.grid()
plt.legend()
plt.show()