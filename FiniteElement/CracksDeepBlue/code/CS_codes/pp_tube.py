#*****************************************************************
#Copyright (c) 2016 Regents of the University of Michigan
#Author: Yan Liu and Matt Collette
#Contact: mdcoll@umich.edu
#See LICENSE files for details of license and code 
#*****************************************************************

import nsga2_michigan_vfo as nsga2
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




post_proc = pp.NSGA_PostProcess('cant_tube_opt.db')
post_proc1 = pp.NSGA_PostProcess('cant_tube_mm.db')
opt_result = post_proc.getFront( 100, 0, [1,2])
opt_result1 = post_proc1.getFront( 100, 0, [1,2])


x = post_proc.getIndVariables(opt_result)
x = np.array(x)
x1 = post_proc.getIndVariables(opt_result1)
x1 = np.array(x1)
sita1 = x[:,0]
sita2 = x[:,1]
sita11 = x1[:,0]
sita21 = x1[:,1]


pf = np.array(opt_result)
pf[:,1]=1/pf[:,1]
pf1 = np.array(opt_result1)
pf1[:,1]=1/pf1[:,1]

pf[:,2] = pf[:,2]*np.pi*0.12
pf1[:,2] = pf1[:,2]*np.pi*0.12


fig = plt.figure()
plt.plot(pf[:,1],pf[:,2],  'r^', mfc='none', label = 'Exact FORM analysis', markersize=11)
plt.plot(pf1[:,1],pf1[:,2],  'ro', mfc='none', label = 'Surrogate method', markersize=11)
plt.ticklabel_format(useOffset=False)
plt.xlabel('$\Delta_1*\Delta_2$', fontsize=20)
plt.ylabel('Volume($cm^3$)', fontsize=20)
plt.legend(loc=4)
plt.xlim([0,100])
plt.show()


