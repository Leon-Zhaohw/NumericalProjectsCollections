#*****************************************************************
#Copyright (c) 2016 Regents of the University of Michigan
#Author: Yan Liu and Matt Collette
#Contact: mdcoll@umich.edu
#See LICENSE files for details of license and code 
#*****************************************************************

import nsga2_michigan_vfo as nsga2
import latinhypercube as latin
import math
import numpy as np
import copy
from scipy.optimize import *
import cProfile
import kriging_hess as kriging
import scipy.stats as stats
import latinhypercube as latin
from operator import itemgetter
import scipy.spatial.distance as dis
from pyre import *
import pareto_front_comparison as comp
import post_process as pp 
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt



post_proc = pp.NSGA_PostProcess('CFRP_opt.db')
post_proc1 = pp.NSGA_PostProcess('CFRP_mm.db')
opt_result = post_proc.getFront( 100, 0, [1,2])
opt_result1 = post_proc1.getFront( 100, 0, [1,2])


x = post_proc.getIndVariables(opt_result)
x = np.array(x)
x1 = post_proc.getIndVariables(opt_result1)
x1 = np.array(x1)
sita1 = x[:,7]
sita2 = x[:,8]
sita11 = x1[:,7]
sita21 = x1[:,8]


pf = np.array(opt_result)
pf[:,1]=1/pf[:,1]
pf1 = np.array(opt_result1)
pf1[:,1]=1/pf1[:,1]

pf[:,2] = pf[:,2]*0.001
pf1[:,2] = pf1[:,2]*0.001


fig = plt.figure()
plt.plot(pf[:,1],pf[:,2],  'r^', mfc='none', label = 'Exact FORM analysis', markersize=11)
plt.plot(pf1[:,1],pf1[:,2],  'ro', mfc='none', label = 'Surrogate method', markersize=11)
plt.ticklabel_format(useOffset=False)
plt.xlabel('normalized $\Delta_1*\Delta_2$', fontsize=20)
plt.ylabel('mass($kg$)', fontsize=20)
plt.legend(loc=4)
#plt.xlim([0,100])
plt.show()


fig = plt.figure()
ax = fig.add_subplot(111, projection = '3d')
ax.scatter(sita1, sita2, pf[:,2], c='b', marker = '^')
ax.set_xlabel(r'$ \Delta_1$', fontsize=20)
ax.set_ylabel(r'$ \Delta_2$', fontsize=20)
ax.set_zlabel('mass($kg$)', fontsize=20)
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
plt.ticklabel_format(useOffset=False)
#plt.zticks(fontsize = 20)
plt.ylim([0,0.1])
plt.show()