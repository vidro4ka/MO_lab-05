import random
import sys
import matplotlib.pyplot as plt
import numpy as np
import math

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

    def printer(self, opt_lambda, omega, delta, J):
        if self.radios == 3:
            print("+", 7 * '-', "+", 10 * '-', "+", 22 * "-", "+", 9 * '-', "+", 9 * '-', "+", sep="")
            print("|\th\t|\tdis    |\t\talpha\t\t  |\tw\t\t|\td\t  |")
            print("+", 7 * '-', "+", 10 * '-', "+", 22 * "-", "+", 9 * '-', "+", 9 * '-', "+", sep="")
            print(f"|  {self.lambdas[0]}", end="")
            print(f"  |{round(self.dist[0], 5):>9f}", end="")
            print(f" | {self.optimal_a[0][0]:>0.4f} {self.optimal_a[0][1]:>0.4f} {self.optimal_a[0][2]:>0.4f} ",
                  end="")
            print(f"|{self.omegas[0]:>0.5f} |", end="")
            print(f" {self.deltas[0]:>0.5f} |")
            for i in range(1, len(self.lambdas)):
                print(f"|  {self.lambdas[i]}", end="")
                print(f"  |{round(self.dist[i], 5):>9f}", end="")
                print(f" | {self.optimal_a[i][0]:>0.4f} {self.optimal_a[i][1]:>0.4f} {self.optimal_a[i][2]:>0.4f} ",
                      end="")
                print(f"| {self.omegas[i]:>0.5f} |", end="")
                print(f" {self.deltas[i]:>0.5f} |")
            print("+", 7 * '-', "+", 10 * '-', "+", 22 * "-", "+", 9 * '-', "+", 9 * '-', "+", sep="")

        if self.radios == 5:
            print("+", 7 * '-', "+", 10 * '-', "+", 22 * "-", "+", 9 * '-', "+", 9 * '-', "+", sep="")
            print("|\th\t|\tdis    |\t\talpha\t\t  |\tw\t\t|\td\t  |")
            print("+", 7 * '-', "+", 10 * '-', "+", 22 * "-", "+", 9 * '-', "+", 9 * '-', "+", sep="")
            print(f"|  {self.lambdas[0]}", end="")
            print(f"  |{round(self.dist[0], 5):>9f}", end="")
            print(f" | {self.optimal_a[0][0]:>0.4f} {self.optimal_a[0][1]:>0.4f} {self.optimal_a[0][2]:>0.4f} ",
                  end="")
            print(f"| {self.omegas[0]:>0.5f} |", end="")
            print(f" {self.deltas[0]:>0.5f} |")
            for i in range(1, len(self.lambdas)):
                print(f"|  {self.lambdas[i]}", end="")
                print(f"  |{round(self.dist[i], 5):>9f}", end="")
                print(f" | {self.optimal_a[i][0]:>0.4f} {self.optimal_a[i][1]:>0.4f} {self.optimal_a[i][2]:>0.4f} ",
                      end="")
                print(f"| {self.omegas[i]:>0.5f} |", end="")
                print(f" {self.deltas[i]:>0.5f} |")
            print("+", 7 * '-', "+", 10 * '-', "+", 22 * "-", "+", 9 * '-', "+", 9 * '-', "+", sep="")

    def scater_plot(self):
        temp_color = ["red", "orange", "green", "blue", "pink", "yellow", "lime", "grey", "black", "cyan", "olive"]
        for i in range(0, 11):
            plt.scatter(0, 0, color="purple")
            plt.scatter(self.omegas[i], self.deltas[i], color=temp_color[i])

        plt.legend([
            f'h   |  delta   |  omega',
            f'{self.lambdas[0]}|{self.deltas[0]:>0.5f}|{self.omegas[0]:>0.5f}',
            f'{self.lambdas[1]}|{self.deltas[1]:>0.5f}|{self.omegas[1]:>0.5f}',
            f'{self.lambdas[2]}|{self.deltas[2]:>0.5f}|{self.omegas[2]:>0.5f}',
            f'{self.lambdas[3]}|{self.deltas[3]:>0.5f}|{self.omegas[3]:>0.5f}',
            f'{self.lambdas[4]}|{self.deltas[4]:>0.5f}|{self.omegas[4]:>0.5f}',
            f'{self.lambdas[5]}|{self.deltas[5]:>0.5f}|{self.omegas[5]:>0.5f}',
            f'{self.lambdas[6]}|{self.deltas[6]:>0.5f}|{self.omegas[6]:>0.5f}',
            f'{self.lambdas[7]}|{self.deltas[7]:>0.5f}|{self.omegas[7]:>0.5f}',
            f'{self.lambdas[8]}|{self.deltas[8]:>0.5f}|{self.omegas[8]:>0.5f}',
            f'{self.lambdas[9]}|{self.deltas[9]:>0.5f}|{self.omegas[9]:>0.5f}',
            f'{self.lambdas[10]}|{self.deltas[10]:>0.5f}|{self.omegas[10]:>0.5f}',
        ])
        print("+-------+---------+---------+---------+")
        print("|   h*  |    J    |    w    |     d   |")
        print("+-------+---------+---------+---------+")
        print(f"|  {optimal_lambda:>0.1f}", end="")
        print(f"  | {optimal_J:>0.5f}", end="")
        print(f" | {optimal_omega:>0.5f} |", end="")
        print(f" {optimal_delta:>0.5f} |")
        print("+-------+---------+---------+---------+")
        print()
        plt.xlabel("omega")
        plt.ylabel("delta")
        plt.grid()
        plt.show()


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
    ob_3.printer(optimal_lambda, optimal_omega, optimal_delta, optimal_J)
    filter_k = np.array([float((k * math.pi) / 100) for k in range(1, 100)])
    plt.plot(k, main_signal, k, main_noise, filter_k, filtering)
    plt.legend(['f(x) = sin(x) + 0.5', 'Noise', 'Filtering'])
    plt.xlabel('x')
    plt.ylabel('f(x)')
    plt.title('Sliding Window size r = 3')
    plt.grid()
    plt.show()
    ob_3.scater_plot()

    r = 5
    ob_5 = SolverClass(main_signal, main_noise, lambda_l, r)
    ob_5.method(main_noise, lambda_l)
    ob_5.distance()
    optimal_J, optimal_alpha, optimal_lambda, optimal_delta, optimal_omega, min_dist, filtering = ob_5.filter()
    ob_5.printer(optimal_lambda, optimal_omega, optimal_delta, optimal_J)
    filter_k = np.array([float((k * math.pi) / 100) for k in range(1, 98)])
    plt.plot(k, main_signal, k, main_noise, filter_k, filtering)
    plt.legend(['f(x) = sin(x) + 0.5', 'Noise', 'Filtering'])
    plt.xlabel('x')
    plt.ylabel('f(x)')
    plt.title('Sliding Window size r = 5')
    plt.grid()
    plt.show()
    ob_5.scater_plot()
