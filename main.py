import random
import sys

import matplotlib.pyplot as plt
import numpy as np
import math


def plot_func(xk, y, colour):
    plt.plot(xk, y, colour)
    plt.title("functions")
    plt.xlabel("x")
    plt.ylabel("f(x)")
    plt.grid(True, alpha=1, linestyle='--')


def gen_sigma_and_noise(a, K, y):
    sigma_k = np.array(float(random.uniform(-a, a)))
    for i in range(1, K):
        sigma_k = np.append(sigma_k, float(random.uniform(-a, a)))
    fnk = np.array(y[0] + sigma_k[0])
    for i in range(1, K):
        fnk = np.append(fnk, y[i] + sigma_k[i])
    return fnk


def gen_alpha(r, M):
    alfa_arr = []

    for i in range(r):
        alfa_arr.append(0)
    alfa_arr[M] = round(random.uniform(0, 1), 4)

    for m in range(M, r - 2):
        summary = 0
        for s in range(m, r - m):
            summary += alfa_arr[s]
        alfa_arr[m - 1] = alfa_arr[r - m] = round(0.5 * random.uniform(0, 1 - summary), 4)

    summary = 0
    for s in range(1, r - 1):
        summary += alfa_arr[s]
    alfa_arr[0] = alfa_arr[r - 1] = round(0.5 * (1 - summary), 4)

    return alfa_arr


def arithmetic_mean(M, noise_fnk, alfa_a):
    func_mov_av = []
    for k in range(M, 100 - M + 1):
        summary = 0
        for j in range(k - M, k + M + 1):
            summary += noise_fnk[j - 1] * alfa_a[j + M + 1 - k - 1]
        func_mov_av.append(summary)
    return func_mov_av


def omegasearcher(func_mov_av, M):
    summary = 0
    for k in range(1, 101 - 2 * M):
        summary += abs(func_mov_av[k] - func_mov_av[k - 1])
    return summary


def deltasearcher(func_mov_av, func_noise, M):
    summary = 0
    for k in range(0, 101 - 2 * M):
        summary += abs(func_mov_av[k] - func_noise[k])
    summary = summary / 100
    return summary


def dist(func_mov_av, func_noise):
    return abs(omegasearcher(func_mov_av)) + abs(deltasearcher(func_mov_av, func_noise))


def signal(xk):
    return np.sin(xk) + 0.5


class SolverClass(object):

    def __init__(self, _signal, _noise, _lambda, _radius):
        self.optimal_a = []
        self.omegas = []
        self.deltas = []
        self.J_min = []
        self.dist = []
        self.radios = _radius
        self.signal = _signal
        self.noise = _noise
        self.lambdas = _lambda

    def method(self, noise_array, lambda_array):
        M = int((self.radios - 1) / 2)
        N = int((math.log(1 - 0.95)) / math.log(1 - (0.01 / math.pi)))
        for i in range(len(lambda_array)):
            best_alpha = np.array([])
            best_omega = 0
            best_delta = 0
            J_min = sys.maxsize
            for V in range(0, N):
                alpha_array = gen_alpha(self.radios, M)
                average_moving = arithmetic_mean(M, noise_array, alpha_array)
                omega = omegasearcher(average_moving, M)
                delta = deltasearcher(average_moving, self.noise, M)
                J = lambda_array[i] * omega + (1 - lambda_array[i]) * delta
                if J < J_min:
                    J_min = J
                    best_alpha = alpha_array.copy()
                    best_omega = omega
                    best_delta = delta
                average_moving.clear()
                self.J_min.append(J_min)
                self.optimal_a.append(best_alpha)
                self.omegas.append(best_omega)
                self.deltas.append(best_delta)

    def distance(self):
        for i in range(0, len(self.deltas)):
            self.dist.append(abs(self.omegas[i]) + self.deltas[i])

    def filter(self):
        M = int((self.radios - 1) / 2)
        min_dist = min(self.dist)
        index = self.dist.index(min_dist)
        optimal_omega = self.omegas[index]
        optimal_delta = self.deltas[index]
        optimal_alpha = self.optimal_a[index]
        optimal_lambda = 0.1 * index
        optimal_J = optimal_lambda * optimal_omega + (1 - optimal_lambda) * optimal_delta
        filtering = []
        for k in range(M, 100 - M + 1):
            summary = 0
            for j in range(k - M, k + M + 1):
                summary += main_noise[j - 1] * optimal_alpha[j + M + 1 - k - 1]
            filtering.append(summary)
        return optimal_J, optimal_alpha, optimal_lambda, optimal_delta, optimal_omega, min_dist, filtering


if __name__ == "__main__":
    k = np.array([float((k * math.pi) / 100) for k in range(0, 101)])
    main_signal = np.array([signal(elem) for elem in k])
    main_noise = np.array([signal(elem) + (random.random() / 2 - 0.25) for elem in k])
    lambda_l = np.array([l_elem / 10 for l_elem in range(0, 11)])

    r = 3
    ob_3 = SolverClass(main_signal, main_noise, lambda_l, r)
    ob_3.method(main_noise, lambda_l)
    ob_3.distance()
    optimal_J, optimal_alpha, optimal_lambda, optimal_delta, optimal_omega, min_dist, filtering = ob_3.filter()
    filter_k = np.array([float((k * math.pi) / 100) for k in range(1, 100)])
    plt.plot(k, main_signal, k, main_noise, filter_k, filtering)
    plt.legend(['f(x) = sin(x) + 0.5', 'Noise', 'Filtering'])
    plt.xlabel('x')
    plt.ylabel('f(x)')
    plt.title('Размер скользящего окна r = 3')
    plt.grid()
    plt.show()

    r = 5
    ob_3 = SolverClass(main_signal, main_noise, lambda_l, r)
    ob_3.method(main_noise, lambda_l)
    ob_3.distance()
    optimal_J, optimal_alpha, optimal_lambda, optimal_delta, optimal_omega, min_dist, filtering = ob_3.filter()
    filter_k = np.array([float((k * math.pi) / 100) for k in range(1, 98)])
    plt.plot(k, main_signal, k, main_noise, filter_k, filtering)
    plt.legend(['f(x) = sin(x) + 0.5', 'Noise', 'Filtering'])
    plt.xlabel('x')
    plt.ylabel('f(x)')
    plt.title('Размер скользящего окна r = 5')
    plt.grid()
    plt.show()