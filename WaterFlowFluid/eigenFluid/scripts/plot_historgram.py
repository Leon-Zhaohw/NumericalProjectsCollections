import numpy as np
import matplotlib.pyplot as plt
from matplotlib import pyplot as pl
import matplotlib.pyplot as plt

def read_from_file(fname):
    error = []
    f = open(fname)
    for line in f:
        temp = line.split()
        error.append((float)(temp[0]))
    f.close()
    error = np.asarray(error)
    return error


eig1600 = read_from_file("./hist_T2_8162.txt")
eig100 = read_from_file("./hist_T2_3577.txt")
eig0 = read_from_file("./hist_T2_1010.txt")

#convert to percentage.
eig1600 = eig1600 * 100.0
eig100 = eig100 * 100.0
eig0 = eig0 * 100.0

#error = error*100.0
#error1 = error1*100.0

index = np.arange(0, eig1600.size)

font_size = 30
linge_width = 2
#rint time_total_rnd.shape
#print edge_discover_rnd.shape

fig, ax = plt.subplots()
#ax.set_ylabel('Singular values ', fontsize = font_size)
#ax.set_xlabel('', fontsize = font_size)
#ax.set_title('Singular values for different tensor slice', fontsize = font_size )

plt.xlim(0, 20)
#ax.set_title('number of edges found vs time: Tower of London', fontsize = font_size)
#print(error.shape)
#ax.plot( color = 'blue', lw = linge_width)
#ax.plot( color = 'green', lw = linge_width)
#ax.plot(time_total_rnd.transpose(), edge_discover_rnd.transpose(), color = 'purple', lw = linge_width)
#ax.plot(time_total_bsl.transpose(), edge_discover_bsl.transpose(), color = 'red', lw = linge_width)

line1, = ax.plot(index, eig1600.transpose(), color = 'blue' ,label="n = 8162", lw = linge_width)
line2, = ax.plot(index, eig100.transpose(),color = 'green' ,label = "n = 3577", lw = linge_width)
line3, = ax.plot(index, eig0.transpose(),color = 'red' ,label = "n = 1010", lw = linge_width)

#ax.legend(loc='upper left', bbox_to_anchor=(0., 1.0), fontsize = font_size)

#plt.xticks(np.arange(0, x_ceil, 200))
for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(font_size)
for tick in ax.yaxis.get_major_ticks():
               tick.label.set_fontsize(font_size)   
plt.show()