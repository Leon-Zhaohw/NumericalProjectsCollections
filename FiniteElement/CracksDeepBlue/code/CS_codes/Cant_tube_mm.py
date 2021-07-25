# -*- coding: utf-8 -*-
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
    
    return Sy - sigma_max
    
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
    




class surrogate:
    def __init__(self, function, samplesize, lb, ub):
        self.func = function
        self.size = samplesize
        self.lb = lb
        self.ub = ub
        self.var = len(lb)
        self.x,self.y = self.sampling()
        self.sita = np.matrix(np.ones((1,self.var)))
        self.model = kriging.Kriging(self.x,self.y,self.sita)
        self.model.solve()

        
    def sampling(self):
        dist = stats.uniform
        #dist = stats.beta
        pars = (0,1)
        #pars = (1,5) #beta
        samp=latin.lhs(dist, pars,(self.size,self.var))
        for i in range(self.var):
            time = self.ub[i] - self.lb[i]
            samp[:,i] = samp[:,i]*time+self.lb[i]
        y = np.zeros(self.size).reshape(self.size,1)
        for i in range(self.size):
            y[i] = self.func(samp[i])
        return samp, y 
        
        
    def add(self, xnew, ynew):
        
        self.x = np.vstack((self.x, xnew))
        self.y = np.vstack((self.y, ynew))
        self.model = kriging.Kriging(self.x,self.y,self.sita)
        self.model.solve()
       
    def predict(self, xpred):
        ypred , ypred_j, ypred_h, mse = self.model.predictor(np.matrix(xpred))
        return ypred, ypred_j, ypred_h, mse
        
        




class Test_Problem_mm(nsga2.Problem):

    
    def __init__(self, offset=0):
        genelowbound = [0.0, 5.0, 3.0, 38.0]
        geneupbound = [10.0, 15.0, 6.0,  44.0]
        nsga2.Problem.__init__(self,2,0,4,genelowbound, geneupbound)     
        self._offset = offset

        self.sg = surrogate(reliability, 150,  \
                                      [0.0, 5.0, 3.0, 38.0], \
                                      [10.0, 15.0, 6.0,  44.0])
          
        self.count = 0
        self.store = []

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
    
    
    def kriging_w(self, x, sita1, sita2):
        [s1,s2, t, d] = x
        ypred, ypred_j, ypred_h , mse = self.sg.predict(x)
        J_I = ypred_j[0, 0:2]
        H_I =  ypred_h[0:2, 0:2]
        a =   J_I*np.linalg.inv(H_I)
        I_w = [s1 - a[0,0],  s2 - a[0,1]]
        if I_w[0]>sita1:
            I_w[0] = sita1
        elif I_w[0]<=0.0:
            I_w[0] = 0.0
        if I_w[1]>sita2:
            I_w[1] = sita2
        elif I_w[1]<=5.0:
            I_w[1] = 5.0
        x_w = [I_w[0], I_w[1], t, d]
        return x_w , ypred , mse     
     
    def constraint(self, ind, metamodel=None):
        t = ind.chromosome[2]                  
        d = ind.chromosome[3]
        sita1 = ind.chromosome[0]
        sita2 = ind.chromosome[1]
        beta_t = 3.0
        self.count += 1
        x_w , Y_0, s0 = self.kriging_w([sita1, sita2, t,d], sita1, sita2)
        x_w1, Y_1, s1 = self.kriging_w(x_w, sita1, sita2)
        Y_2, y2_j, y2_h, s2 = self.sg.predict(x_w1)
        if (abs(Y_2[0,0]-Y_1[0,0])>= 1e-3):
            #print Y_2, Y_1
            #print self.count
            self.store.append(self.count)
            y_w = reliability(x_w1)
            self.sg.add(x_w1, y_w)
            x_w , Y_w, m1= self.kriging_w(x_w1, sita1, sita2)
            index_w, y_w_j, y_w_h, error = self.sg.predict(x_w)
        else:
            x_w = x_w1; index_w = Y_2; error = s2
        

        if error[0,0] == 0.0:
            p = max(beta_t - index_w[0,0], 0)
        else:            
            u = (index_w[0,0] - beta_t)/(error[0,0])
        
            if (index_w[0,0]>=beta_t) & (u>=3.0):
                p = 0.0
            elif (index_w[0,0]>=beta_t) & (u<3.0):
                print 'add', index_w
                y_a = reliability(x_w)
                self.sg.add(x_w, y_a)
                p = max(0.0, beta_t - y_a)
            else:
                p = beta_t - index_w[0,0]
        return([p], None)


def run_mm():
    myProb = Test_Problem_mm()
    myOpt = nsga2.Optimizer(myProb, 0.8, 4.0, 0.01, 4.0)
    myOpt.run('cant_tube_mm.db', 'three_obj_test', 1007, 100, 100) 

run_mm()
#cProfile.run('run_mm()')

#def test():
#    start_seed = 1001
#    numseeds = 10
#    pop = 100
#    gen = 100 
#    for seedoffset in range(0, numseeds):
#        seed = start_seed + seedoffset
#        current_problem = Test_Problem_mm()
#        filename = 'tube_mmm' + str(seed) + '.db'
#        opt = nsga2.Optimizer(current_problem, 0.8, 4.0, 0.01, 4.0)
#        opt.run(filename, "Tube " + 
#                str(seed), seed, gen, pop)
#        print('Finished run ' + filename)
#test()

