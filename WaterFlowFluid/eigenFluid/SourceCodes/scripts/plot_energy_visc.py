import numpy as np
import matplotlib.pyplot as plt
from matplotlib import pyplot as pl
import matplotlib.pyplot as plt

def read_from_file(fname):
    e = []
    f = open(fname)
    for line in f:
        temp = line.split()
        e.append((float)(temp[0]))
    f.close()
    e = np.asarray(e)
    #ed = np.asarray(ed)
    result = {}
    result['e'] = e
    #result['ed'] = ed
    return result

data0 = read_from_file("./E_visc0.txt")['e']
data0 = data0[2:-1]
data0 = data0 / data0[0] * 100.0

data1 = read_from_file("./E_visc0.0001.txt")['e']
data1 = data1[2:-1]
data1 = data1 / data1[0] * 100.0

data2 = read_from_file("./E_visc0.0002.txt")['e']
data2 = data2[2:-1]
data2 = data2 / data2[0] * 100.0

index = np.arange(0, len(data1))

font_size = 30
linge_width = 2

fig, ax = plt.subplots()
#ax.set_ylabel('%  error ', fontsize = font_size)
#ax.set_xlabel('Timestep', fontsize = font_size)
#ax.set_title('Error comparison, visc = 0.002', fontsize = font_size )
num = len(index)

line1, = ax.plot(index, data1.transpose(), color = 'blue' , lw = linge_width)
line2, = ax.plot(index, data2.transpose(), color = 'green' , lw = linge_width)
line3, = ax.plot(index, data0.transpose(), color = 'red' , lw = linge_width)
#line4, = ax.plot(index, derivative.transpose()*0.05, color = 'red', lw = linge_width)
#line3, = ax.plot([1,600],[0,0], color = 'black')
plt.ylim(0, 110)
ax.legend(loc='upper left', bbox_to_anchor=(0., 1.0), fontsize = font_size)

#plt.xticks(np.arange(0, x_ceil, 200))
for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(font_size)
for tick in ax.yaxis.get_major_ticks():
               tick.label.set_fontsize(font_size)   
plt.show()