import numpy as np
import matplotlib.pyplot as plt
from matplotlib import pyplot as pl
import matplotlib.pyplot as plt

def read_from_file(fname):
    w = []
    num = []
    f = open(fname)
    for line in f:
        temp = line.split()
        w.append((float)(temp[0]))
        num.append((float)(temp[1]))
    f.close()
    w = np.asarray(w)
    num = np.asarray(num)
    result = {}
    result['w'] = w
    result['num'] = num
    return result

data = read_from_file("../weight.txt")

index = np.arange(0, len(data['w']))

font_size = 30
linge_width = 2

fig, ax = plt.subplots()
ax.set_ylabel('Basis weight ', fontsize = font_size)
ax.set_xlabel('Wavenumber |K|^2', fontsize = font_size)
ax.set_title('Forward scattering amplified 5.0', fontsize = font_size )
num = len(index)

line1, = ax.plot(data['w'].transpose(), data['num'].transpose(), color = 'blue' , lw = linge_width)
#plt.ylim(-1.2,0)
#ax.legend(loc='upper left', bbox_to_anchor=(0., 1.0), fontsize = font_size)

#plt.xticks(np.arange(0, x_ceil, 200))
for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(font_size)
for tick in ax.yaxis.get_major_ticks():
               tick.label.set_fontsize(font_size)   
plt.show()