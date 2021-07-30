import numpy as np
import matplotlib.pyplot as plt
from matplotlib import pyplot as pl
import matplotlib.pyplot as plt

def read_from_file(fname):
    el = []
    ed = []
    f = open(fname)
    for line in f:
        temp = line.split()
        el.append((float)(temp[0]))
        ed.append((float)(temp[1]))
    f.close()
    el = np.asarray(el)
    ed = np.asarray(ed)
    result = {}
    result['el'] = el
    result['ed'] = ed
    return result

data = read_from_file("../test.txt")

index = np.arange(0, len(data['el']))

font_size = 30
linge_width = 2

fig, ax = plt.subplots()
#ax.set_ylabel('%  error ', fontsize = font_size)
ax.set_xlabel('steps', fontsize = font_size)
#ax.set_title('Error comparison, visc = 0.002', fontsize = font_size )
num = len(index)
derivative = np.zeros(data['el'].shape)
dx = 1.0 / num
for i in xrange(0, num - 1):
    derivative[i] = (data['el'][i+1] - data['el'][i])/ dx

line1, = ax.plot(index, data['el'].transpose(), color = 'blue' , lw = linge_width)
line2, = ax.plot(index, data['ed'].transpose()*2, color = 'green' , lw = linge_width)
line4, = ax.plot(index, derivative.transpose()*0.05, color = 'red', lw = linge_width)
line3, = ax.plot([1,600],[0,0], color = 'black')
#plt.ylim(-0.01,0.03)
ax.legend(loc='upper left', bbox_to_anchor=(0., 1.0), fontsize = font_size)

#plt.xticks(np.arange(0, x_ceil, 200))
for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(font_size)
for tick in ax.yaxis.get_major_ticks():
               tick.label.set_fontsize(font_size)   
plt.show()