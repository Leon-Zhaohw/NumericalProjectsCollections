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


error = read_from_file("./coefficients_ziperr_old.txt")
error1 = read_from_file("./coefficients_zip1err_old.txt")

error = error*100.0
error1 = error1*100.0

index = np.arange(0, 398)

font_size = 30
linge_width = 2
#rint time_total_rnd.shape
#print edge_discover_rnd.shape

fig, ax = plt.subplots()
ax.set_ylabel('%  error ', fontsize = font_size)
ax.set_xlabel('steps', fontsize = font_size)
ax.set_title('Error comparison, visc = 0.002', fontsize = font_size )

#ax.set_title('number of edges found vs time: Tower of London', fontsize = font_size)
print(error.shape)
#ax.plot( color = 'blue', lw = linge_width)
#ax.plot( color = 'green', lw = linge_width)
#ax.plot(time_total_rnd.transpose(), edge_discover_rnd.transpose(), color = 'purple', lw = linge_width)
#ax.plot(time_total_bsl.transpose(), edge_discover_bsl.transpose(), color = 'red', lw = linge_width)

line1, = ax.plot(index, error.transpose(), color = 'blue' ,label="%75 compression", lw = linge_width)
line2, = ax.plot(index, error1.transpose(),color = 'green' ,label = "%50 compression", lw = linge_width)
ax.legend(loc='upper left', bbox_to_anchor=(0., 1.0), fontsize = font_size)

#plt.xticks(np.arange(0, x_ceil, 200))
for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(font_size)
for tick in ax.yaxis.get_major_ticks():
               tick.label.set_fontsize(font_size)   
plt.show()