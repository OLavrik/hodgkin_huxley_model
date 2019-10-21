import flanagan.integration.DerivnFunction;
import flanagan.integration.RungeKutta;
import org.math.plot.Plot2DPanel;

import javax.swing.*;

public class Model {
    final static double g_l_max = 0.003; //  ms/mm^2
    final static double g_k_max = 0.36;
    final static double g_na_max = 1.2;
    final static double e_l_max = -54.387;  // mV
    final static double e_k_max = -77;
    final static double e_na_max = 50;
    final static double dt = 0.01; // ms
    double min_delta_t; // ms
    double max_delta_t;
    double min_i;
    double max_i;
    double V = -64.9964;
    int t;
    double[] I;
    int[] t_I;
    int size;
    double[] V_time;//=new double[(int) (this.size/this.dt)];
    double[] I_time;
    double[] t_time;

    double m = 0.0530;
    double h = 0.5960;
    double n = 0.3177;

    double g_na;
    double g_k;
    double[] m_time;
    double[] h_time;
    double[] n_time;

    int num = 0;

    class Dervi implements DerivnFunction {

        @Override
        public double[] derivn(double t, double[] y) {
            double[] dydt = new double[4];

            dydt[0] = I_time[num] - g_na_max * Math.pow(m, 3.0) * h * (y[0] - e_na_max) - g_k_max * Math.pow(n, 4.0) * (y[0] - e_k_max) - g_l_max * (y[0] - e_l_max); //V
            dydt[1] = 0.1 * (y[0] + 40) / (1 - Math.exp(-(y[0] + 40) / 10)) * (1.0 - y[1]) - 4 * Math.exp(-0.0556 * (y[0] + 65)) * y[1]; //m
            dydt[2] = 0.07 * Math.exp(-0.05 * (y[0] + 65)) * (1.0 - y[2]) - 1 / (1 + Math.exp(-0.1 * (y[0] + 35))) * y[2]; //h
            dydt[3] = 0.01 * (y[0] + 55) / (1 - Math.exp(-(y[0] + 55) / 10)) * (1.0 - y[3]) - 0.125 * Math.exp(-(y[0] + 65) / 80) * y[3]; //n

            return dydt;

        }
    }


    public void plot_graph(double[] x, double[] y) {

        Plot2DPanel plot = new Plot2DPanel();


        // add a line plot to the PlotPanel
        plot.addLinePlot("my plot", x, y);

        // put the PlotPanel in a JFrame like a JPanel
        JFrame frame = new JFrame("a plot panel");
        frame.setSize(600, 600);
        frame.setContentPane(plot);
        frame.setVisible(true);

    }


    public void plot() {
        this.plot_graph(this.t_time, this.V_time);
        plot_graph(t_time, h_time);
        plot_graph(t_time, n_time);
        plot_graph(t_time, m_time);
        plot_graph(t_time, I_time);
        plot_graph(V_time, m_time);
        plot_graph(V_time, n_time);
        plot_graph(V_time, h_time);
    }

    public void calculat_with_dif() {
        V_time[0] = V;
        m_time[0] = m;
        h_time[0] = h;
        n_time[0] = n;
        for (int i = 0; i < size - 1; i++) {
            dif(i);
        }
    }


    public void calculat() {

        for (int i = 0; i < size; i++) {

            g_na = g_na_max * (Math.pow(m, 3.0)) * h;
            g_k = g_k_max * (Math.pow(n, 4.0));
            double g_sum = g_k + g_na + g_l_max;
            double v_inf = ((g_na * e_na_max + g_k * e_k_max + g_l_max * e_l_max) + I_time[i]) / g_sum;
            double tau_v = .01 / g_sum;
            V = v_inf + (V - v_inf) * Math.exp(-dt / tau_v);
            double alpha_m = 0.1 * (V + 40) / (1 - Math.exp(-(V + 40) / 10));
            double beta_m = 4 * Math.exp(-0.0556 * (V + 65));
            double alpha_n = 0.01 * (V + 55) / (1 - Math.exp(-(V + 55) / 10));
            double beta_n = 0.125 * Math.exp(-(V + 65) / 80);
            double alpha_h = 0.07 * Math.exp(-0.05 * (V + 65));
            double beta_h = 1 / (1 + Math.exp(-0.1 * (V + 35)));
            double tau_m = 1 / (alpha_m + beta_m);
            double tau_h = 1 / (alpha_h + beta_h);
            double tau_n = 1 / (alpha_n + beta_n);
            double m_inf = alpha_m * tau_m;
            double h_inf = alpha_h * tau_h;
            double n_inf = alpha_n * tau_n;
            this.m = m_inf + (this.m - m_inf) * Math.exp(-dt / tau_m);
            this.h = h_inf + (this.h - h_inf) * Math.exp(-dt / tau_h);
            this.n = n_inf + (this.n - n_inf) * Math.exp(-dt / tau_n);
            V_time[i] = V;
            m_time[i] = m;
            h_time[i] = h;
            n_time[i] = n;
            t_time[i] = i;

        }
    }

    public void make_I_array() {
        int end = 0;
        this.I_time = new double[(int) ((this.t + 100 * 4) / this.dt) + 200];
        for (int i = 0; i < this.I.length; i++) {
            for (int j = 0; j < (int) (this.t_I[i] / this.dt); j++) {

                this.I_time[end + j] = this.I[i];
            }
            end = (int) ((this.t_I[i] + this.t) / this.dt) + end;
        }


    }


    public void generate_args(int size) {
        this.I = new double[size];
        this.t_I = new int[size];
        for (int i = 0; i < size; i++) {
            I[i] = Math.random() * ((this.max_i - this.min_i) + 1) + this.min_i;
            t_I[i] = (int) (Math.random() * ((this.max_delta_t - this.min_delta_t) + 1) + this.min_delta_t);

        }


    }


    Model() {
        this.min_delta_t = 20.0;  // ms
        this.max_delta_t = 100.0;
        this.min_i = 5.0 * 0.001;
        this.max_i = 100.0 * 0.001;
        this.t = 50;
        generate_args(4);
        make_I_array();
        this.size = (int) ((this.t + 100 * 4) / this.dt);

        this.m_time = new double[size];
        this.h_time = new double[size];
        this.n_time = new double[size];
        this.V_time = new double[size];
        this.t_time = new double[size];
    }


    Model(int mode2) {
        this.min_delta_t = 1.0;  // ms
        this.max_delta_t = 1.0;
        this.t = 1;
        this.min_i = 50.0;
        this.max_i = 200.0;
        generate_args(2);
        make_I_array();
        this.size = (int) ((this.t + 100 * 4) / this.dt);

        this.m_time = new double[size];
        this.h_time = new double[size];
        this.n_time = new double[size];
        this.V_time = new double[size];
        this.t_time = new double[size];
    }

    public void dif(int i) {
        Dervi dn = new Dervi();
        int size_array = 4;
        double hi = 0.01;
        double x0 = (double) (i);              // initial value of x
        double xn = (double) (i) + 1;
        double[] y0 = new double[size_array];         // initial values of the y[i]
        double[] yn = new double[size_array];
        y0[0] = V_time[i];
        y0[1] = m_time[i];
        y0[2] = h_time[i];
        y0[3] = n_time[i];
        RungeKutta rk = new RungeKutta();
        rk.setInitialValueOfX(x0);
        rk.setFinalValueOfX(xn);
        rk.setInitialValuesOfY(y0);
        rk.setStepSize(hi);
        yn = rk.fourthOrder(dn);

        rk.setToleranceScalingFactor(1e-10);
        rk.setToleranceAdditionFactor(1e-8);


        this.V_time[i + 1] = yn[0];
        this.m_time[i + 1] = yn[0];
        this.h_time[i + 1] = yn[0];
        this.n_time[i + 1] = yn[0];

        System.out.println("Fourth order Runge-Kutta procedure");
        System.out.println("The value of y[0] at x = " + xn + " is " + yn[0]);
        System.out.println("The value of y[1] at x = " + xn + " is " + yn[1]);
        System.out.println("The value of y[2] at x = " + xn + " is " + yn[2]);
        System.out.println("The value of y[2] at x = " + xn + " is " + yn[3]);
        System.out.println("Number of iterations = " + rk.getNumberOfIterations());


    }

}


