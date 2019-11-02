import random as rd
import numpy as np
from membrane import Membrane
import matplotlib.pyplot as plt

# import pylab as plt
color = ['g', 'r', 'b']


class Model:
    def plot(self,  y_arr,x_arr, ylabel, xlabel,fun=1):
        print("meow")
        plt.figure()
        if (fun>1):
            for i in range(y_arr.shape[0]):
                plt.plot( x_arr,y_arr[i], color[i])
        else:
            plt.plot(x_arr,y_arr,'k')
        plt.ylabel(ylabel)
        plt.xlabel(xlabel)
        plt.show()

    def run(self):
        m = Membrane(self.size, self.I_time, self.I_part)
        # self.plot(m.V, m.time, "V(mV)", "t(mc)")
        print(np.stack((m.m, m.n, m.h)).shape)
        i_inj = [m.i_run(i) for i in m.time]
        self.plot( np.stack((m.m, m.n, m.h)),m.time, "m/n/h", "t(mc)",3)
        self.plot(i_inj, m.time, "n(mV)", "t(mc)")
        self.plot(m.V, m.time, "h(mV)", "t(mc)")

    def generate_args(self, count):
        self.size = 0
        for i in range(count):
            self.I_part.append(round(rd.random() * ((self.max_i - self.min_i) + 1) + self.min_i, 2))
            self.I_time.append([round(self.size, 2), round(
                self.size + rd.random() * ((self.max_delta_t - self.min_delta_t) + 1) + self.min_delta_t, 2)])
            self.size = round(self.I_time[i][1] + self.delta_t, 2)
        self.size = self.size - self.delta_t
        print(self.size)

    def __init__(self, mode=1):
        self.I_time = []
        self.I_part = []
        if mode == 1:
            self.min_delta_t = 20.0  # ms
            self.max_delta_t = 100.0;
            self.min_i = 5.0;
            self.max_i = 100.0;
            self.delta_t = 50;
            self.generate_args(4);
        else:
            self.min_delta_t = 1.0  # ms
            self.max_delta_t = 1.0;
            self.min_i = 50.0;
            self.max_i = 200.0;
            self.delta_t = 1;
            self.generate_args(2);


if __name__ == '__main__':
    m = Model()
    m.run()
