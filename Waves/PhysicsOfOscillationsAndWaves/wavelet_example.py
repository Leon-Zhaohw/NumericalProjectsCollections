'''
Example of how to do wavelet analysis in python
NB NB NB This thing either gets stuck or does not work!!

The scritp was in 2017 translated by Sebastian G. Winther-Larsen from a matlab script originally written by Arnt Inge Vistnes
'''

from matplotlib import pyplot as plt
import numpy as np 
import pywt
import sounddevice as sd
import soundfile as sf

# Scalogram function
def scalogram(data):
    bottom = 0

    vmin = np.min(map(lambda x: np.min(abs(x)), data))
    vmax = np.max(map(lambda x: np.min(abs(x)), data))

    plt.gca().set_autoscale_on(False)

    for row in range(len(data)):
        scale = 2.0 ** (row - len(data))

        plt.imshow(
            np.array([abs(data[row])]),
            interpolation = 'nearest',
            vmin = vmin,
            vmax = vmax,
            extent = [0 , 1, bottom, bottom + scale]
        )

        bottom += scale

# Load signal
data, Fs = sf.read("svarttrost2.wav")
tree = pywt.wavedec(data, 'db5')

# Plotting
scalogram(data)
plt.show()