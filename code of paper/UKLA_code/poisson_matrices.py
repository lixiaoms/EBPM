import numpy as np
import scipy.linalg as scl
import finitediff_minus_one as fmo


def omega(Z):
    return np.exp(Z)


def omega_inv(X):
    return np.log(X)


# Variational model
def shrink(Y, lbd=1, T=100, centering='both'):
    m = Y.shape[0]
    k = Y.shape[1]

    X_MLplus = np.maximum(Y, 0.5)

    # Rule of thumb: lambda should be approximatively proportional
    # to the standard deviation of Y which is sqrt(Y)
    lbd = lbd * np.sqrt(Y.mean())

    t = 1
    gamma = 1 / Y.max()
    Z = omega_inv(X_MLplus)
    for i in range(T):
        Z_old = Z.copy()

        grad_F_Z = omega(Z) - Y
        Z = Z - gamma * grad_F_Z
        if lbd > 0:
            if centering is 'row':
                mean_Z = Z.mean(axis=0, keepdims=True)
            if centering is 'col':
                mean_Z = Z.mean(axis=1, keepdims=True)
            if centering is 'both':
                mean_Z = \
                    Z.mean(axis=0, keepdims=True) + \
                    Z.mean(axis=1, keepdims=True) - \
                    Z.mean()
            if centering is None:
                mean_Z = 0
            Z = Z - mean_Z
            U, S, VH = np.linalg.svd(Z, full_matrices=False)
            S = np.maximum(S - gamma * lbd, 0)
            Z = np.dot(U * S, VH)
            Z = mean_Z + Z

        # FISTA update
        t_old = t
        t = (1 + np.sqrt(1 + 4 * t**2)) / 2.
        Z = Z + (t_old - 1) / t * (Z - Z_old)

    X_hat = omega(Z)
    return X_hat

def complete(Y, obs, lbd=1, T=100, centering='both'):
    m = Y.shape[0]
    k = Y.shape[1]

    X_MLplus = np.maximum(Y, 0.5)

    # Rule of thumb: lambda should be approximatively proportional
    # to the standard deviation of Y which is sqrt(Y)
    lbd = lbd * np.sqrt(np.sum(Y*obs)/np.sum(obs))

    t = 1
    gamma = 1 / Y.max()
    Z = omega_inv(X_MLplus)
    for i in range(T):
        Z_old = Z.copy()

        grad_F_Z = np.multiply(omega(Z) - Y, obs)
        Z = Z - gamma * grad_F_Z
        if lbd > 0:
            if centering is 'row':
                mean_Z = Z.mean(axis=0, keepdims=True)
            if centering is 'col':
                mean_Z = Z.mean(axis=1, keepdims=True)
            if centering is 'both':
                mean_Z = \
                    Z.mean(axis=0, keepdims=True) + \
                    Z.mean(axis=1, keepdims=True) - \
                    Z.mean()
            if centering is None:
                mean_Z = 0
            Z = Z - mean_Z
            U, S, VH = np.linalg.svd(Z, full_matrices=False)
            S = np.maximum(S - gamma * lbd, 0)
            Z = np.dot(U * S, VH)
            Z = mean_Z + Z

        # FISTA update
        t_old = t
        t = (1 + np.sqrt(1 + 4 * t**2)) / 2.
        Z = Z + (t_old - 1) / t * (Z - Z_old)

    X_hat = omega(Z)
    return X_hat


def KL(X, X_hat_func, Y):
    X_hat = X_hat_func(Y)
    return (X_hat - X * np.log(X_hat)).sum()


def UKL(X_hat_func, Y, eps=.1, directions=None, order=1):
    if not directions:
        directions = []
        for o in range(order):
            directions = directions + \
                [(2 * np.random.binomial(1, .5, size=y.shape) - 1)]
    X_hat = X_hat_func(Y)
    log_X_hat_m1 = fmo.findiff_m1_expansion(lambda Y: np.log(X_hat_func(Y)),
                                            Y, eps, directions[:order+1],
                                            order=order+1)
    return (X_hat - Y * log_X_hat_m1).sum()
