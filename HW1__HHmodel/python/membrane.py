import random
import numpy as np
from scipy.integrate import odeint


class Membrane:
    g_l_max = 0.3  # ms / mm ^ 2
    g_k_max = 36
    g_na_max = 120
    e_l_max = -54.387  # mV
    e_k_max = -77
    e_na_max = 50
    dt = 0.01  # ms
    c_m=0.01

    def alpha_m(self, voltage):
        return 0.1 * (voltage + 40.0) / (1.0 - np.exp(-(voltage + 40.0) / 10.0))

    def beta_m(self, V):
        return 4.0 * np.exp(-(V + 65.0) / 18.0)

    def alpha_h(self, V):
        return 0.07 * np.exp(-(V + 65.0) / 20.0)

    def beta_h(self, V):
        return 1.0 / (1.0 + np.exp(-(V + 35.0) / 10.0))

    def alpha_n(self, V):
        return 0.01 * (V + 55.0) / (1.0 - np.exp(-(V + 55.0) / 10.0))

    def beta_n(self, V):
        return 0.125 * np.exp(-(V + 65) / 80.0)

    def i_na(self, V, m, h):
        return self.g_na_max * m ** 3 * h * (V - self.e_na_max)

    def i_k(self, V, n):
        return self.g_k_max * n ** 4 * (V - self.e_k_max)

    def i_l(self, V):
        return self.g_l_max * (V - self.e_l_max)

    def i_run(self, t):
        for i in range(len(self.I_time)):
            if self.I_time[i][0] <= t <= self.I_time[i][1]:
                return self.I_part[i]
        return 0

    @staticmethod
    def dif(X, t, self):

        V, m, h, n = X

        dVdt = (self.i_run(t) - self.i_na(V, m, h) - self.i_k(V, n) - self.i_l(V)) / self.c_m
        dmdt = self.alpha_m(V) * (1.0 - m) - self.beta_m(V) * m
        dhdt = self.alpha_h(V) * (1.0 - h) - self.beta_h(V) * h
        dndt = self.alpha_n(V) * (1.0 - n) - self.beta_n(V) * n
        return dVdt, dmdt, dhdt, dndt

    def solve(self):

        x0 = [-65, 0.05, 0.6, 0.32]
        x = odeint(self.dif, x0, self.time, args=(self,))
        self.V = x[:, 0]
        self.m = x[:, 1]
        self.h = x[:, 2]
        self.n = x[:, 3]


    def __init__(self, size, I_time, I_part):
        self.time = np.arange(0, size, 0.01)
        self.I_time = I_time
        self.I_part = I_part
        self.solve()
