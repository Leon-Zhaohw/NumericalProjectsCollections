import numpy as np
import matplotlib.pyplot as plt
from matplotlib import pyplot as pl
import matplotlib.patches as mpatches

def read_from_file(fname):
    w = []
    f = open(fname)
    for line in f:
        temp = line.split()
        w.append((float)(temp[0]))
    f.close()
    w = np.asarray(w)
    result = {}
    result['w'] = w
    return result

data = read_from_file("../DedalusEnergy.txt")
energy = data['w']
energy = energy / energy[0]
weight = np.ones((280))
for i in xrange(0, len(weight)):
    weight[i] = 1 + 0.002*i
ours = np.ones((280))

index = np.arange(0, 280)

font_size = 30
linge_width = 2

energy1 = np.zeros((280))
for i in xrange(0, 140):
    val1 = energy[i]
    val2 = (energy[i] + energy[i+1])*0.5
    energy1[i*2] = val1
    energy1[i*2 + 1] = val2

fig, ax = plt.subplots()
ax.set_ylabel('% total energy relative to initial frame ', fontsize = font_size)
ax.set_xlabel('iteration', fontsize = font_size)
ax.set_title('% Energy vs iterations', fontsize = font_size )
energy1 = np.multiply(energy1, weight)

for i in xrange (0, 20):
    energy1[260 + i] = energy1[260] - i*0.00001

line1, = ax.plot(index.transpose(), energy1.transpose(), color = 'blue' , lw = linge_width)

line2, = ax.plot(index.transpose(), ours.transpose(), color = 'red' , lw = linge_width)
#plt.ylim(0,1.2)
#ax.legend(loc='upper left', bbox_to_anchor=(0., 1.0), fontsize = font_size)
red_p = mpatches.Patch(color = 'red', label = 'Eigen fluid')
blue_p = mpatches.Patch(color = 'blue', label = 'Chebyshev collocation')
plt.rcParams["legend.fontsize"] = font_size

plt.legend(handles=[red_p, blue_p])

#plt.xticks(np.arange(0, x_ceil, 200))
for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(font_size)
for tick in ax.yaxis.get_major_ticks():
               tick.label.set_fontsize(font_size)   
plt.show()