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
import post_process as pp 
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from pyre import *


__author__ = "Yan Liu, Han Koo Jeong and Matthew Collette"
__copyright__ = "Copyright 2015, The Regents of the University of Michigan"
__license__ = "See license.dat in source distribution"
__status__ = "Development"



def weight(x):
    [thk_w,  thk_cr, thk_p,  Tr_thk, ths_cr, ths_w, Ts_thk] = x
    #thk_p = 3*thk_w
    #thk_cr = thk_w
    density = 1600.0     #! CFRP mass density (kg/m^3)
    den_form= 20.0       #! Former mass density (kg/m^3)
    #load_P  = 0.0364      #! Applied load (MPa)
    len_l   = 4200.0     #! Length of a plated grillage (mm)
    len_b   = 3200.0     #! Width of a plated grillage (mm)    
    #thk_cr = thk_w
#! It is assumed that b_r = c_r = b_s = c_s = thickness (thk)        
    #Tr_thk  = 180.0       #! Nominal depth of longitudinal beams      
    #PI = 3.1415926
    #E       = 140000.0   #! CFRP Elastic modulus (MPa)
    density = 1600.0     #! CFRP mass density (kg/m^3)
    den_form= 20.0       #! Former mass density (kg/m^3)
    #load_P  = 0.0364      #! Applied load (MPa)
    len_l   = 4200.0     #! Length of a plated grillage (mm)
    len_b   = 3200.0     #! Width of a plated grillage (mm)        
#! It is assumed that b_r = c_r = b_s = c_s = thickness (thk)        
    Ts_thk  = 180.0       #! Nominal depth of longitudinal beams 
    #       #! Nominal depth of transverse beams    
 #! Input variables -> -> -> -> -> -> -> -> -> -> -> -> -> -> -> -> -> -> -> -> -> -> ->
    r       = 3          #! Number of longitudinal beams
    s       = 5          #! Number of transverse beams
    X       = len_l*2/4 #! Coordinate value alongside longitudinal beams (mm); mid-point
    Y       = len_b*3/6  #! Coordinate value alongside transverse beams (mm); mid-point
    spac_r  = 800.0      #! Longitudinal beam frame spacing
    spac_s  = 700.0      #! Transverse beam frame spacing
    Pth     = 2          #! The Pth longitudinal stiffener
    Qth     = 3          #! The Qth transverse stiffener
    
#! It is assumed that the effective breadth of longitudinal and transverse beams 
#! is equal to 20% of their beam spacings.
#! It is assumed that the width of flange is equal to 30% of its effective breadth  
#! It is assumed that the thickness of plating-flange is equal to three times that of 
#! 'thk' for top-hat stiffened single skin panel.   
    b_r    = thk_cr              #! Thickness of the flange of longitudinal beam (mm)
    c_r    = thk_w              #! Thickness of the web of longitudinal beam	(mm)
    f_r    = thk_p            #! Thickness of the plating-flange of longitudinal beam (mm)
    d_r    = Tr_thk-(b_r+f_r) #! Length of the web of longitudinal beam (mm)
    e_r    = spac_r*0.2       #! Effective breadth of longitudinal beam (mm)
    a_r    = e_r*0.3          #! Width of the flange of longitudinal beam (mm)   
    b_s    = ths_cr              #! Thickness of the flange of transverse beam (mm)
    c_s    = ths_w           #! Thickness of the web of transverse beam (mm)
    f_s    = f_r            #! Thickness of the plating-flange of transverse beam (mm)
    d_s    = Ts_thk-(b_s+f_s) #! Length of the web of transverse beam (mm)
    e_s    = spac_s*0.2       #! Effective breadth of transverse beam (mm)
    a_s    = e_s*0.3          #! Width of the flange of transverse beam (mm)
    W_panel  = len_l*len_b*thk_p*density*0.000001   #! Panel weight (g)
    W_F_lon  = (a_r*b_r+2*(c_r*d_r))*len_l*density*0.000001*r   #! Weight of the longitudinal stiffeners (g)
    W_F_tra  = (a_s*b_s+2*(c_s*d_s))*len_b*density*0.000001*s   #! Weight of the transverse stiffeners (g)
    W_Former_lon	= (a_r-2*c_r)*d_r*len_l*den_form*0.000001*r   #! Weight of the former of the longitudinal stiffeners (g)
    W_Former_tra = (a_s-2*c_s)*d_s*len_b*den_form*0.000001*s   #! Weight of the former of the transverse stiffeners (g)
    W_Frame_Overlap = (r*s)*(((d_r*2*c_r)+(a_r*b_r))*a_s*density*0.000001+   \
        ((a_r-2*c_r)*d_r*a_s*den_form*0.000001))   #! Weight of overlapped top-hat stiffeners (g)

    W_Frame  = W_F_lon+W_F_tra+W_Former_lon+W_Former_tra-W_Frame_Overlap   #! Total weight of the stiffeners	(g)
    W_Total  = W_panel+W_Frame  # ! Total weight of the plated grillage (g) 
    #Str_Msth_Ten = M_sth*Yc_s/I_s
    return W_Total


def limitstate(X1,X2,X3, X4, X5, X6, X7, X8, X9):  
    #[E, load_P, thk_w, thk_p, thk_cr,  D_total] = x
    # Input is Thickness    
    #thk = 25.0  
    # X1 = Elastic modulus    X2 = Applied load
    # X3 , X4,  X5, X6 = thk_w, thk_p, thk_cr,  D_total     
    # X7 : computation model uncertainty
    E = X1; load_P = X2; thk_w = X3
    thk_cr = X4; thk_p = X5; ths_w = X6
    ths_cr = X7
    Tr_thk = X8     
    Ts_thk = X9

    #thk_cr = thk_w
    #thk_p = thk_w*3
    ''' 
    INTEGER :: M, N, r, s, Pth, Qth, N_panels
    REAL, PARAMETER :: PI = 3.141592
    REAL    :: E, load_P, len_l, len_b, X, Y
    REAL    :: spac_r, spac_s, Tr_thk, Ts_thk, thk
    REAL    :: a_r, b_r, c_r, d_r, e_r, f_r, Yc_r, I_r
    REAL    :: a_s, b_s, c_s, d_s, e_s, f_s, Yc_s, I_s
    REAL    :: Coeff_W(20,20), p_W, W
    REAL    :: Coeff_M_rth(20,20), p_M_rth, M_rth
    REAL    :: Coeff_M_sth(20,20), p_M_sth, M_sth
    REAL    :: Str_Mrth_Ten, Str_Mrth_Com
    REAL    :: Str_Msth_Ten, Str_Msth_Com
    REAL    :: density, W_panel, W_F_lon, W_F_tra
    REAL    :: den_form, W_Former_lon, W_Former_tra
    REAL    :: W_Frame, W_Frame_Overlap, W_Total    
    '''
    #Tr_thk = 180.0
    
    
    PI = 3.1415926
    #E       = 140000.0   #! CFRP Elastic modulus (MPa)
    density = 1600.0     #! CFRP mass density (kg/m^3)
    den_form= 20.0       #! Former mass density (kg/m^3)
    #load_P  = 0.0364      #! Applied load (MPa)
    len_l   = 4200.0     #! Length of a plated grillage (mm)
    len_b   = 3200.0     #! Width of a plated grillage (mm)
        
#! It is assumed that b_r = c_r = b_s = c_s = thickness (thk)    
    
    #Tr_thk  = 100.0       #! Nominal depth of longitudinal beams 
    #Ts_thk  = Tr_thk       #! Nominal depth of transverse beams
    
 #! Input variables -> -> -> -> -> -> -> -> -> -> -> -> -> -> -> -> -> -> -> -> -> -> ->

    r       = 3          #! Number of longitudinal beams
    s       = 5          #! Number of transverse beams
    X       = len_l*2/4 #! Coordinate value alongside longitudinal beams (mm); mid-point
    Y       = len_b*3/6  #! Coordinate value alongside transverse beams (mm); mid-point
    spac_r  = 800.0      #! Longitudinal beam frame spacing
    spac_s  = 700.0      #! Transverse beam frame spacing
    Pth     = 2          #! The Pth longitudinal stiffener
    Qth     = 3          #! The Qth transverse stiffener
    
#! It is assumed that the effective breadth of longitudinal and transverse beams 
#! is equal to 20% of their beam spacings.
#! It is assumed that the width of flange is equal to 30% of its effective breadth  
#! It is assumed that the thickness of plating-flange is equal to three times that of 
#! 'thk' for top-hat stiffened single skin panel.
   
    b_r    = thk_cr              #! Thickness of the flange of longitudinal beam (mm)
    c_r    = thk_w              #! Thickness of the web of longitudinal beam	(mm)
    f_r    = thk_p            #! Thickness of the plating-flange of longitudinal beam (mm)
    d_r    = Tr_thk-(b_r+f_r) #! Length of the web of longitudinal beam (mm)
    e_r    = spac_r*0.2       #! Effective breadth of longitudinal beam (mm)
    a_r    = e_r*0.3          #! Width of the flange of longitudinal beam (mm)
   
    b_s    = ths_cr             #! Thickness of the flange of transverse beam (mm)
    c_s    = ths_w             #! Thickness of the web of transverse beam (mm)
    f_s    = f_r            #! Thickness of the plating-flange of transverse beam (mm)
    d_s    = Ts_thk-(b_s+f_s) #! Length of the web of transverse beam (mm)
    e_s    = spac_s*0.2       #! Effective breadth of transverse beam (mm)
    a_s    = e_s*0.3          #! Width of the flange of transverse beam (mm)
		   
    Yc_r = (a_r*b_r*(b_r/2)+2*(c_r*d_r*(b_r+d_r/2))+ \
           e_r*f_r*(b_r+d_r+f_r/2))/(a_r*b_r+2*c_r*d_r+e_r*f_r)      
           	 #! The centroid of the area for the longitudinal beams
   
    Yc_s = (a_s*b_s*(b_s/2)+2*(c_s*d_s*(b_s+d_s/2))+e_s*f_s*(b_s+d_s+f_s/2))/(a_s*b_s+2*c_s*d_s+e_s*f_s)     
    	 #! The centroid of the area for the transverse beams
   
    I_r = (a_r*b_r**3/12)+a_r*b_r*(Yc_r-b_r/2)**2+           \
      2*(c_r*d_r**3/12+c_r*d_r*(b_r+d_r/2-Yc_r)**2)+       \
        (e_r*f_r**3/12)+e_r*f_r*(b_r+d_r+f_r/2-Yc_r)**2   
 
#! The second moment of the area for the longitudinal beams
   
    I_s = (a_s*b_s**3/12)+a_s*b_s*(Yc_s-b_s/2)**2+           \
      2*(c_s*d_s**3/12+c_s*d_s*(b_s+d_s/2-Yc_s)**2)+       \
        (e_s*f_s**3/12)+e_s*f_s*(b_s+d_s+f_s/2-Yc_s)**2   
 
#! The second moment of the area for the transverse beams
   
    W = 0
    M_rth = 0
    M_sth = 0
   
    m = np.arange(1,16,2)
    
    n = np.arange(1,16,2)
    for i in range(8):
        for j in range(8):
            M = m[i]
            N = n[j]
            
            CoeffW_M_N = (16*load_P*len_l**4*len_b/(E*I_r))/(PI**6*M*N*(M**4*(r+1)+(I_s/I_r)*(len_l**3/len_b**3)*N**4*(s+1)))      
                                      
            p_W = CoeffW_M_N*np.sin((M*PI*X/len_l))*np.sin((N*PI*Y/len_b))
            W = W + p_W	       #! Deflection of the plated grillage

            Coeff_M_rth_M_N = (16*load_P*len_l**2*len_b)/((PI**4*N/M)*(M**4*(r+1)+(I_s/I_r)*(len_l**3/len_b**3)*N**4*(s+1)))                       
                                
            p_M_rth = Coeff_M_rth_M_N*np.sin((M*PI*X/len_l))*np.sin((N*PI*Pth/(r+1)))
            M_rth = M_rth + p_M_rth	  # ! Bending moment for the Pth longitudinal stiffener
   
            Coeff_M_sth_M_N = (16*load_P*len_b**2*len_l)/((PI**4*M/N)*(N**4*(s+1)+(I_r/I_s)*(len_b**3/len_l**3)*M**4*(r+1)))                      
                      
            p_M_sth = Coeff_M_sth_M_N*np.sin((M*PI*Qth/(s+1)))*np.sin((N*PI*Y/len_b))
            M_sth = M_sth + p_M_sth     #! Bending moment for the Qth transverse stiffener


    
    return 4200.0/150 - W #composite_test.Grillage_Structures_Tophat_Singleskin(X2,X3, X4, X5)*X6
 



   
def reliability(x):
    [ thk_w,thk_cr, thk_p,ths_w, ths_cr, Tr_thk, Ts_thk, COV1, COV2] = [x[0],x[1],x[2], x[3], x[4], x[5], x[6], x[7], x[8]]
    limit_state = LimitState(limitstate)
    options = AnalysisOptions()
    options.printResults(False)
    stochastic_model = StochasticModel()
    # Define random variables
    stochastic_model.addVariable( Normal('X1',140000, COV1*140000.0) )
    stochastic_model.addVariable( Normal('X2',0.0364,  COV2*0.0364) )      
    stochastic_model.addVariable( Normal('X3',thk_w,0.03*thk_w) )
    stochastic_model.addVariable( Normal('X4',thk_cr,0.03*thk_cr) )
    stochastic_model.addVariable( Normal('X5',thk_p,0.03*thk_p) )  
    stochastic_model.addVariable( Normal('X6',ths_w,0.03*ths_w) )
    stochastic_model.addVariable( Normal('X7',ths_cr,0.03*ths_cr) )
    stochastic_model.addVariable( Normal('X8',Tr_thk,0.03*Tr_thk) )
    stochastic_model.addVariable( Normal('X9',Ts_thk,0.03*Ts_thk) )
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
        pars = (0,1)
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






class Test_Problem(nsga2.Problem):
    
    def __init__(self, offset=0):
        #thk_w,thk_cr, thk_p,ths_w, ths_cr, Tr_thk, Ts_thk, COV1, COV2
        genelowbound = [1.0,1.8, 5.0, 1.0, 1.8, 150, 150, 0.010001,  0.010001]
        geneupbound = [15.0, 15.0, 45.0, 15.0, 15.0,  180, 180, 0.1, 0.1]
        nsga2.Problem.__init__(self,2,0,9,genelowbound, geneupbound)  
        self.mm = surrogate(reliability, 200,  [1.0,1.8, 5.0, 1.0, 1.8, 150, 150,  0.05,  0.05],  \
                                          [15.0, 15.0, 45.0, 15.0, 15.0,  180, 180, 0.25,  0.25])        

        self._offset = offset

    def function1(self, chromosome):

        delta1 = chromosome[7]
        delta2 = chromosome[8]
        a1 = (delta1-0.01)/(0.1-0.01)
        a2 = (delta2-0.01)/(0.1-0.01)
        return 1/(a1*a2)
            
 
    def function2(self,chromosome):
        thk_w = chromosome[0]
        thk_cr = chromosome[1]
        thk_p = chromosome[2]
        ths_w = chromosome[3]
        ths_cr = chromosome[4]
        Tr_thk = chromosome[5]
        Ts_thk = chromosome[6]
        return weight([thk_w,thk_cr, thk_p, ths_w, ths_cr, Tr_thk, Ts_thk])
               
    def Eval(self,ind, metamodel=None):       
        f1 = self.function1(ind.chromosome) 
        f2 = self.function2(ind.chromosome)
        f = [f1,f2]
        return (f, None)
    
    def kriging_w(self, x, sita):
        [a1,a2,a3,a4,a5,a6,a7, F1, F2] = x
        [sita1, sita2] = sita
        ypred, ypred_j, ypred_h , mse = self.mm.predict(x)
        J_I = ypred_j[0, 7:9]
        H_I =  ypred_h[7:9, 7:9]
        a =   J_I*np.linalg.inv(H_I)
        I_w = [F1 - a[0,0],  F2 - a[0,1]]
        if I_w[0]>0.15+sita1:
            I_w[0] = 0.15+sita1
        elif I_w[0]<=0.15-sita1:
            I_w[0] = 0.15-sita1
        if I_w[1]>0.15+sita2:
            I_w[1] = 0.15+sita2
        elif I_w[1]<=0.15-sita2:
            I_w[1] =0.15-sita2
        x_w = [I_w[0], I_w[1]]
        return x_w , ypred , mse      
        
    def constraint(self, ind, metamodel=None):
        thk_w =  ind.chromosome[0]
        thk_cr = ind.chromosome[1]
        thk_p =  ind.chromosome[2]
        ths_w =  ind.chromosome[3]
        ths_cr = ind.chromosome[4]
        Tr_thk = ind.chromosome[5]
        Ts_thk = ind.chromosome[6]
        delta1 = ind.chromosome[7]
        delta2 = ind.chromosome[8]
        beta_w = 3.0
        x_w, ypred, mse = self.kriging_w([thk_w,thk_cr, thk_p,ths_w, ths_cr, Tr_thk, Ts_thk,0.15,0.15],[delta1,delta2])
        x_w1, ypred1, mse1 = self.kriging_w([thk_w,thk_cr, thk_p,ths_w, ths_cr, Tr_thk, Ts_thk,x_w[0],x_w[1]],[delta1,delta2]) 
        x_w2, ypred2, mse2 = self.kriging_w([thk_w,thk_cr, thk_p,ths_w, ths_cr, Tr_thk, Ts_thk,x_w1[0],x_w1[1]],[delta1,delta2])  
           
        
        if ypred2[0,0]+mse2[0,0] <= beta_w:
            p1 = beta_w - ypred2[0,0]
             
        else:
            #print x_w, x_w1
            diff = abs(x_w1[0]-x_w2[0])+abs(x_w1[1]-x_w2[1])  
            if diff >=1e-3:           
                y_a = reliability([thk_w,thk_cr, thk_p,ths_w, ths_cr, Tr_thk, Ts_thk,x_w2[0],x_w2[1]])
                self.mm.add([thk_w,thk_cr, thk_p,ths_w, ths_cr, Tr_thk, Ts_thk,x_w2[0],x_w2[1]],y_a)
                x_w, ypred, mse = self.kriging_w([thk_w,thk_cr, thk_p,ths_w, ths_cr, Tr_thk, Ts_thk,x_w2[0],x_w2[1]],[delta1,delta2])
            else:
                x_w = x_w2
                
            ypred_w, ypred_jw, ypred_hw , mse1= self.mm.predict([thk_w,thk_cr, thk_p,ths_w, ths_cr, Tr_thk, Ts_thk,x_w[0],x_w[1]])    
            index_w = ypred_w[0,0]
        
            if mse1[0,0]==0.0:
                U = 3.50
            else:
                U = (index_w-beta_w)/mse1[0,0]
            
            if (index_w>= beta_w) and U>=3.0 :
                p1 = 0.0
            elif (index_w>=beta_w) and U<3.0:
                #print 'update'
                y_a = reliability([thk_w,thk_cr, thk_p,ths_w, ths_cr, Tr_thk, Ts_thk,x_w2[0],x_w2[1]])
                self.mm.add([thk_w,thk_cr, thk_p,ths_w, ths_cr, Tr_thk, Ts_thk,x_w2[0],x_w2[1]], y_a)
                index_w = y_a
                p1 = max(0.0, beta_w-y_a)
        
            else:
                p1 = beta_w - index_w
        
        return([p1], None)
        
        
        
def run():
    myProb = Test_Problem()
    myOpt = nsga2.Optimizer(myProb, 0.8, 4.0, 0.01, 4.0)
    myOpt.run('CFRP_mm.db', 'CFRP_test', 1024, 100, 100)  
    
#cProfile.run('run()')

run()
