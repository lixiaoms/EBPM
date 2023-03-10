import numpy as np
from functools import reduce


# Estimate directional derivatives of desired orders
# with desired directions at x
def findiff_dir_der(f, x, eps, directions, order=1):
    if order < 0:
        raise NameError('order cannot be negative')

    if len(directions) != order:
        raise NameError('the number of directions must be equal the order')

    if order == 0:
        return f(x)

    diff = findiff_dir_der(f, x + eps * directions[0],
                           eps, directions[1:], order=order-1) - \
        findiff_dir_der(f, x - eps * directions[0],
                        eps, directions[1:], order=order-1)

    return diff / (2 * eps)


# Estimate f(x - e_i)_i for all i based on Taylor expansion with desired orders
# on desired directions and finite differences with step eps
def findiff_m1_expansion(f, x, eps, directions, order=1):
    res = 0
    for o in range(order+1):
        eps_order = eps**(1 / np.maximum(o, 1)) / 4
        Jx = findiff_dir_der(f, x, eps_order, directions[:o], order=o)
        prod = reduce(lambda x, y: x * y, directions[:o], 1)
        res += (-1)**o * prod * Jx / np.math.factorial(o)
    return res
