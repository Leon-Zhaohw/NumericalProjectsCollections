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

import post_process as pp 
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt


def limitstate(X1,X2,X3,X4,X5,X6,X7,X8,X9, X10, X11):
    t = X1; d = X2; L1 = X3; L2 = X4; F1 = X5*1000; F2 = X6*1000
    P = X7*1000; T = X8*1000; Sy = X9   
    sita1 = X10*np.pi/180; sita2 = X11*np.pi/180
    A = (np.pi/4)*(d**2 - (d-2*t)**2)
    I = (np.pi/64)*(d**4 - (d-2*t)**4)    
    tor = T*d/(4*I)
    
    M = F1*L1*np.cos(sita1)+F2*L2*np.cos(sita2)
    sigma_x = (P+F1*np.sin(sita1)+F2*np.sin(sita2))/A + M*(d/2)/I
    sigma_max = (sigma_x**2 + 3*tor**2)**0.5
    
    return Sy - sigma_max+ np.random.normal(0.0, 0.001)
    
def reliability(x):
    [sita1, sita2, t,  d] = [x[0],x[1], x[2], x[3]]
    limit_state = LimitState(limitstate )
    options = AnalysisOptions()
    options.printResults(False)
    stochastic_model = StochasticModel()
    # Define random variables
    stochastic_model.addVariable( Normal('X1',t,0.1) )    
    stochastic_model.addVariable( Normal('X2',d,0.5) )   
    stochastic_model.addVariable( Uniform('X3',120,0.25) )
    stochastic_model.addVariable( Uniform('X4',60.0,0.25) )
    stochastic_model.addVariable (Normal('X5',3.0,0.3))
    stochastic_model.addVariable (Normal('X6',3.0,0.3))
    stochastic_model.addVariable (Gumbel('X7',12.0,1.2))
    stochastic_model.addVariable (Normal('X8',90.0,9.0))
    stochastic_model.addVariable (Normal('X9',220,22.0))
    stochastic_model.addVariable (Normal('X10',sita1,0.0001))       
    stochastic_model.addVariable (Normal('X11',sita2,0.0001))       
    
    Analysis = Form(analysis_options=options, stochastic_model=stochastic_model, limit_state=limit_state)
    
    beta = Analysis.getBeta()
    #failure = Analysis.getFailure()
    return beta 
    





class Test_Problem(nsga2.Problem):
    
    def __init__(self, offset=0):
        genelowbound = [0.0, 5.0, 3.0, 38.0]
        geneupbound = [10.0, 15.0, 6.0,  44.0]
        nsga2.Problem.__init__(self,2,0,4,genelowbound, geneupbound)     
        self._offset = offset
        self.count = 0

    def function1(self, chromosome):
        sita1 = chromosome[0]
        sita2 = chromosome[1]-5.0
        return 1/(sita1*sita2+1e-5)
            
 
    def function2(self,chromosome):
        t = chromosome[2]
        d = chromosome[3]
        return t*(d-t)
               
    def Eval(self,ind, metamodel=None):       
        f1 = self.function1(ind.chromosome) 
        f2 = self.function2(ind.chromosome)
        f = [f1,f2]
        return (f, None)
        
        
    def constraint(self, ind, metamodel=None):
        t = ind.chromosome[2]                
        d = ind.chromosome[3]
        sita1 = ind.chromosome[0]
        sita2 = ind.chromosome[1]
      
        def lsf(x):
            self.count += 1
            if x[0]>=sita1:
                p1 = x[0]-sita1
            else:
                p1 = 0.0
            if x[0]<0.0:
                p2 = 0.0-x[0]
            else:
                p2 = 0.0
            if x[1]>=sita2:
                p3 = x[1]-sita2
            else:
                p3 = 0.0
            if x[1]<5.0:
                p4 = 5.0-x[1]
            else:
                p4 = 0.0    
        
            return reliability([x[0],x[1], t, d])+100000.0*sum([p1,p2,p3,p4])
            
        w_t = optimize.fmin_powell(lsf, [0.0, 5.0], args=(), xtol=0.0001, ftol=0.0001, maxiter=None, maxfun=None, full_output=0, disp=0, retall=0, callback=None, direc=None)    
        index_w = reliability([w_t[0], w_t[1], t, d])
        print index_w
        if index_w>= 3.0:
            p = 0.0
        else:
            p = 3.0 - index_w

        return([p], None)
        
        
        
        
def run():
    myProb = Test_Problem()
    myOpt = nsga2.Optimizer(myProb, 0.8, 4.0, 0.01, 4.0)
    myOpt.run('cant_tube_opt.db', 'three_obj_test', 1024, 100, 100)  

run()


