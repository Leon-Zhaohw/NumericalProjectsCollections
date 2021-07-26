'''
This program simulates slit diffraction experiment.
The aim is to find the intensity distribution of light passing
through a single og a double slit.

The scritp was in 2017 translated by Sebastian G. Winther-Larsen from a matlab script originally written by Arnt Inge Vistnes
'''

import sys
import numpy as np
from matplotlib import pyplot as plt
import time
import progressbar

# GENERAL PARAMETERS
n_per_lambda = 4                        # Points per wavelength
screen_distance = 10000 * n_per_lambda  # Distance between slit and screen
screen_width = 6000                     # Screen width for study in wavelengths
N = screen_width * n_per_lambda         # Number of points to calculate
N_half = N / 2

# EXPERIMENT SPECIFIC PARAMETERS
slit_width = 10                         # Width of slit in wavelengths
slit_separation = 80                    # Separation of slits (if two)

M = 200                         # Full 1/e width of Gaussian amplitude profile
R = 1000                        # Radius of curvature in wavefront


# ALLOCATING ARRAYS FOR COMPUTATION
x  = np.linspace(-screen_width/2, \
                    screen_width/2, N)      # Relative position array for plot
x2 = np.linspace(-N, N, 2*N + 1)            # For check of helper functions
x0 = np.zeros((N, 2))                      # Excitation, amplitude and phase
x1 = np.zeros((N, 2))                      # At screen, amplitude and phase
                   # Distances. Reduction factor
                                            # and phase correction.

# DIFFERENT EXPERIMENT CONFIGURATIONS

# Single slit experiment
def single_slit(x0, n_per_lambda, slit_separation, slit_width, N, N_half):
    print("Single slit experiment.")
    m = slit_width * n_per_lambda / 2
    N_half = int(N_half)
    m = int(m)
    x0[N_half - m : N_half + m - 1, 0] = 1.0
    return x0

# Double slit experiment
def double_slit(x0, n_per_lambda, slit_separation, slit_width, N, N_half):
    print("Double slit experiment.")
    m = slit_width * n_per_lambda
    k_left  = (slit_separation / 2 + slit_width / 2) * n_per_lambda
    k_right = (slit_separation / 2 - slit_width / 2) * n_per_lambda
    N_half  = int(N_half)
    k_left  = int(k_left)
    k_right = int(k_right)
    m = int(m)
    x0[N_half - k_left  + 1 : N_half - k_left  + m,     0] = 1.0
    x0[N_half + k_right + 1 : N_half + k_right + m - 1, 0] = 1.0
    return x0

def straight_edge(x0, n_per_lambda, slit_separation, slit_width, N, N_half):
    print("Straight edge diffraction.")
    x0[int(N/4):int(N)] = 1.0
    return x0

def gaussian_excitation(x0, n_per_lambda, slit_separation, slit_width, N, N_half):
    print("Gaussian excitation.")
    width = 200 * n_per_lambda / 2
    dummy = (np.arange(1, N + 1) - N_half) / width
    dummy = dummy * dummy
    x0[:, 0] = np.exp(-dummy)
    # Phase
    R = 1000 # Curvature radius
    y = np.arange(-N_half, N_half)
    R2 = R * R * n_per_lambda * n_per_lambda
    dist = np.sqrt((y * y) + R2)
    fs = dist % n_per_lambda
    x0[:, 1] = fs * (2*np.pi / n_per_lambda)

    return x0

def read_data_from_file(*args):
    print("Reading excitation data from file.")
    f = input("Please enter file name: ")
    x0 = np.loadtxt(f)
    return x0


# DATA GENERATION FUNCTIONS

# Sin, cos, distances, relative phase differences etc.
def generate_relative_position_data(N, screen_distance, screen_width, wavelength):
    y = np.arange(-N, N+1)
    distance_squared = screen_distance * screen_distance
    y2p = y*y + distance_squared
    rnn = np.sqrt(y2p)
    sin_cos = np.zeros((2*N + 1, 2))
    r = np.zeros((2*N +1, 2))
    sin_cos[:, 0] = screen_distance / rnn
    sin_cos[:, 1] = y / rnn
    r[:, 0] = 1 / np.sqrt(rnn)
    fs = rnn % wavelength
    r[:, 1] = fs*(2*np.pi / wavelength)
    return sin_cos, r

# Adding up all contributions for each point on screen
def summation(x, x0, x1, r, N):

    # Progress bar
    bar = progressbar.ProgressBar()
    phasor = np.zeros(x0.shape)

    for i in bar(range(N)):

        rel_position1 = N + 1  - i
        rel_position2 = rel_position1 + N
        amplitude = x0[:, 0] * r[rel_position1:rel_position2, 0]
        phase = x0[:, 1] - r[rel_position1:rel_position2, 1]
        phasor[:, 0] = amplitude * np.cos(phase)
        phasor[:, 1] = amplitude * np.sin(phase)
        phasor_x = np.sum(phasor[:, 0])
        phasor_y = np.sum(phasor[:, 1])
        x1[i - 1, 0] = np.sqrt(phasor_x*phasor_x + phasor_y*phasor_y)
        x1[i - 1, 1] = np.arctan2(phasor_y, phasor_x)

    return x1

# Plotting  diffraction spectrum intensities
def diffraction_plot(x, x0, x1, experiment):
    x1_squared = x1[:, 0] * x1[:, 0] # Intensity
    scaling_factor = np.max(x1_squared) / 8
    plt.plot(x, x0[:, 0]*scaling_factor, '-r')
    plt.plot(x, x1_squared, '-b')
    plt.xlabel("Position on screen")
    plt.ylabel("Relative intensity")
    plt.title("Experiment: " + experiment)
    plt.show()

def write_to_file(x1):
    np.savetxt('output.dat', x1)

# Something else

if __name__ == '__main__':

    # Dictionary of options
    options = {
        1 : single_slit,
        2 : double_slit,
        3 : straight_edge,
        4 : gaussian_excitation,
        5 : read_data_from_file
    }

    # Printing instructions
    if (len(sys.argv) <= 1 or len(sys.argv) > 2):
        print("\nUsage options: \n")
        for key in options:
            print("{} : {}" .format(key, options[key].__name__))
        print("\nExample: 'python {} 1'\n" .format(__file__))
    else:
        config = int(sys.argv[1])

        # Executing chosen option
        x0 = options[config](x0, n_per_lambda, slit_separation, slit_width, N, N_half)

        # Computing sine, cosine, distances, relative phase shifts etc
        sin_cos, r = generate_relative_position_data(N, screen_distance, screen_width, n_per_lambda)

        # Adding contributions for every pint on screen
        x1 = summation(N, x0, x1, r, N)

        # Plotting
        diffraction_plot(x, x0, x1, options[config].__name__)

        # Write to file
        # write_to_file(x1)