'''
This program is meant to illustrate the difference phase- and
group velocity.

The scritp was in 2017 translated by Sebastian G. Winther-Larsen from a matlab script originally written by Arnt Inge Vistnes
'''

import numpy as np
from matplotlib import pyplot as plt
import matplotlib.animation as animation

# Making spatial wave packet
def phase_group_wavepacket(N, x_max, x_lambda, x_sigma):

    x = np.linspace(0, x_max * (N - 1) / N, N)
    x_center = x_max / 8
    x_freq = 1/ x_lambda
    y = np.cos((x - x_center) * 2*np.pi * x_freq)
    convol = np.exp(-((x - x_center) / x_sigma) * ((x - x_center) / x_sigma))
    z = y * convol
    return x, z

# Fourier Transform of
def phase_group_fft(z, N, x_max):

    z_fft = np.fft.fft(z) / N
    A = 2 * np.abs(z_fft)
    theta = np.arctan2(z_fft.real, z_fft.imag)
    x_sample_freq = N/ x_max
    x_freq = np.linspace(0, x_sample_freq * (N - 1) / N, N)
    k = 2 * np.pi * x_freq

    # In order to pick out i_min, i_max
    # plt.plot(A)
    # plt.show()

    return A, theta, k

def phase_group_omega(i_min, i_max, k, dispersion):

    # Creating dictionary int the range of i_min and i_max
    omega = dict.fromkeys(range(i_min, i_max))

    vf_ref = 400
    k_dom = (i_max - i_min)/2.0 + 1

    if (dispersion == -1):
        delta_t = 0.0070
        for i in omega:
            omega[i] = vf_ref * 1.04 * np.sqrt(k[i]/k_dom)

    if (dispersion == 0):
        delta_t = 0.0070
        for i in omega:
            omega[i] = vf_ref * k[i] / k_dom

    if (dispersion == 1):
        delta_t = 0.0070
        for i in omega:
            omega[i] = vf_ref * 1.1 * (k[i] / k_dom)**1.5

    return omega, delta_t

def phase_group_wave(x, t, N, A, phase, k, omega, i_min, i_max):

    z_recon = np.zeros(N)

    # Want one iteration if i_min == i_max
    if i_min == i_max:
        i = i_min
        arg = k[i] * x - omega[i]*t + phase[i]
        z_recon += A[i] * np.cos(arg)

    for i in range(i_min, i_max):
        arg = k[i] * x - omega[i]*t + phase[i]
        z_recon += A[i] * np.cos(arg)

    return z_recon

def animate(frame, data, lines, points, text):

    animate_wave(frame, data, lines)
    animate_phase(frame, data, points)
    animate_text(frame, data, text)

def animate_phase(frame, data, points):

    delta_t = data[1]
    k       = data[5]
    omega   = data[6]
    i_min   = data[7]
    i_max   = data[8]
    x_max   = data[9]

    t = delta_t * frame

    i_mean = int(np.round((i_min + i_max) / 2))

    x0 = x_max / 8.0
    x_avt = x0 + t*omega[i_mean] / k[i_mean]

    points.set_data([x0, x_avt], [0.25, 0.25])

    return points,

def animate_wave(frame, data, lines):

    # Unpacking data
    # data = [x, delta_t, N, A, phase, k, omega, i_min, i_max, x_max]
    x       = data[0]
    delta_t = data[1]
    N       = data[2]
    A       = data[3]
    phase   = data[4]
    k       = data[5]
    omega   = data[6]
    i_min   = data[7]
    i_max   = data[8]
    x_max   = data[9]

    i_mean = int(np.round(i_min + i_max) / 2)
    t = delta_t * frame
    z_recon0 = phase_group_wave(x, t, N, A, phase, k, omega, i_mean, i_mean)
    z_recon  = phase_group_wave(x, t, N, A, phase, k, omega, i_min, i_max)

    lines[0][0].set_data(x, z_recon)
    lines[1][0].set_data(x, 2.5*z_recon0)
    return lines,

def animate_text(frame, data, text):

    delta_t = data[1]
    k       = data[5]
    omega   = data[6]
    i_min   = data[7]
    i_max   = data[8]
    x_max   = data[9]

    t = delta_t * frame

    i_mean = int(np.round((i_min + i_max)/2))

    x0 = x_max / 8.0
    x_avt = x0 + t*omega[i_mean] / k[i_mean]

    text[0].set_position((3.0, 0.8))
    text[0].set_text("Time: {:<5.2f}".format(t))

    text[1].set_position((3.0, 0.65))
    text[1].set_text("x_ref: {:<5.2f}".format(x_avt))

    return text

if __name__ == '__main__':

    # Choosing dispersion, can be
    # -1: Normal dispersion
    #  0: No dispersion
    # +1: Anomalous dispersion
    dispersion = -1

    # Creating a spatial wave packet
    N = 4000
    x_max = 40
    x_lambda = 1
    x_sigma = 2
    x, z = phase_group_wavepacket(N, x_max, x_lambda, x_sigma)
    plt.plot(x,z)
    plt.show()

    # Spatial frequency analysis
    A, phase, k = phase_group_fft(z, N, x_max)

    # Possibility to pick out areas with high amplitude from freqency plot
    i_min = 23
    i_max = 59

    # Finding omega from dispersion relation
    omega, delta_t = phase_group_omega(i_min, i_max, k, dispersion)

    # Figure to plot on
    fig = plt.figure()
    lines = (plt.plot([]), plt.plot([]))
    points, = plt.plot([], 'rx')
    text = [plt.text(0, 0, ""), plt.text(0, 0, "")]
    plt.xlim(0, x_max)
    plt.ylim(-1.04, 1.04)
    if dispersion == -1:
        plt.title("Normal dispersion")
    if dispersion == 0:
        plt.title("No dispersion")
    if dispersion == +1:
        plt.title("Anomalous dispersion")

    # Animation arguments
    data = [x, delta_t, N, A, phase, k, omega, i_min, i_max, x_max]
    ani_args = (data, lines, points, text)

    # Animate!
    ani = animation.FuncAnimation(fig, animate, frames=200, init_func=None,
                    interval=20, fargs=ani_args)

    # Show figure
    plt.show()