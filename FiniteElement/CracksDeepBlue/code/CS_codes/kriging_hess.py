# -*- coding: utf-8 -*-
#*****************************************************************
#Copyright (c) 2016 Regents of the University of Michigan
#Author: Yan Liu and Matt Collette
#Contact: mdcoll@umich.edu
#See LICENSE files for details of license and code 
#*****************************************************************



import numpy as np
from scipy.linalg import inv, qr, cholesky
from scipy.optimize import *

__author__ = "Jiandao Zhu, Matt Collette, Yan Liu"
__copyright__ = "Copyright 2015, The Regents of the University of Michigan"
__license__ = "See license.dat in source distribution"
__status__ = "Development"

class Struct:
    def __init__(self, **kwds):
        self.__dict__.update(kwds)


class Kriging:
    """
    Class using kriging model to fit the provided data set and predict the unknow
    points
    """

    def __init__(self,xpts,ypts,sita0):
        '''
        Imports initial data and sita (correlation) parameter guess 
        
        Parameters (note these require numpy matrix, not array functions)
        ----------
        xpts:  A numpy matrix of points describing the design points where 
               each row is a point in design space.  Any number of columns
        
        ypts:  A nx1 numpy matrix of function responses at the design points
               current code can only handle a single output
               
        site0:  Numpy array of initial guesses for the scaled-data sita (often
                written theta). As many entries as xpts has columns.
       
        Returns
        -------
        None
        '''
        (xpts,ypts)=self._normalize(xpts,ypts)
        npts = xpts.shape[0]
        nvars=xpts.shape[1]
        #structure with data point and intial guess of the correlation function parameters
        Data=Struct(xpts=xpts,ypts=ypts,sita=sita0,npts=npts,nvars=nvars)
        self._Data=Data
        
    def _normalize(self,x,y):
        '''Normalize data'''
        col=x.shape[1]
        row=x.shape[0]
        mx=list()
        sx=list()
        for i in range(col):
            mx.append(np.mean(x[:,i]))
            sx.append(np.std(x[:,i],ddof=1))
        my=np.mean(y)
        sy=np.std(y,ddof=1)
        self._mx=mx
        self._sx=sx
        self._my=my
        self._sy=sy
        mx_mat=np.mat(np.zeros((row,col)))
        sx_mat=np.mat(np.zeros((row,col)))
        my_mat=np.mat(np.zeros((row,1)))
        sy_mat=np.mat(np.zeros((row,1)))
        for i in range(row):
            mx_mat[i,:]=mx
            sx_mat[i,:]=sx
            my_mat[i,0]=my
            sy_mat[i,0]=sy
        x_nom=(x-mx_mat)/sx_mat
        y_nom=(y-my_mat)/sy_mat
        return x_nom,y_nom

    def _varcorr(self,vals,sita,n=1.0):
        '''
        Returns the component of the correlation matrix for a particular
        variable, defined as a sita(|xi-xj|^2).
        Vals should be a 1xn numpy matrix of the variable positions
        Sita is the scale parameter
        n is the exponent value for the correlation function
           1 - Exponential correlation
           2 - Gaussian correlation (default)
        '''
        npts = vals.shape[0]
        
        #Make a matrix of the array of points repeating
        pts_mat = np.matrix(np.ones((npts,npts)))

        #Insert the points as repeated columns
        for i in xrange(npts):
            pts_mat[:,i] = vals

        #Now calculate the transpose
        out = pts_mat - np.transpose(pts_mat)
        out = np.abs(out)
        out = np.power(out,n)
        out *= sita
        return out

    def _corrmatrix(self,vals):
        '''
        Returns the full correlation matrix
        Vals should be matrix of points, each variable in a column
        site should be an array of values of the sita variable
        '''
        #Assert that we have the shapes we want
        assert vals.ndim == 2, \
          "Error in Kriging correlation matrix, variables not 2-D matrix"
        assert self._Data.sita.ndim == 1, \
           "Error in Kriging correlation matrix, sita not passed in as array"
        assert self._Data.sita.shape[0] == vals.shape[1], \
            "Error in Krining correlation matrix, number of sita != num vals"
        
        
        
        #First find the number of variables
        npts = vals.shape[0]
        nvars = vals.shape[1]

        #Make the final output matrix
        out = np.matrix(np.zeros((npts,npts)))

        #Get the contribution from each variable
        for i in xrange(nvars):
            out = out + self._varcorr(vals[:,i],self._Data.sita[i],2)

        out *= -1.0
        out = np.exp(out)
        
        return out

    
    def _objfun(self,sita):
        """
        minimize the function in terms of sita
        """
        #print sita
        self._Data.sita=sita    
        # if get any negative sita,return a large fi value
        for i in range(self._Data.nvars):
            if self._Data.sita[i]<=0:
                fi=1e16
                return fi
        x=np.matrix(np.ones((self._Data.npts,1)))   
        F=np.matrix(np.ones((self._Data.npts,1)))       
        R=np.matrix(np.zeros((self._Data.npts,self._Data.npts)))
        R = self._corrmatrix(self._Data.xpts)
        mu=(10+self._Data.npts)*np.finfo(float).eps
        R=R+np.eye(self._Data.npts)*mu
        #Cholesky factorization
        C=cholesky(R)
        
        # Get least squares solution
        C=np.transpose(C)
        
        Ft=inv(C)*F
        
        (Q,R)=qr(Ft)
        
        Q=Q[:,0]
        
    
        G=np.transpose(Ft)/Q
      
        G=G[0,0]
        Yt=inv(C)*self._Data.ypts
        beta=(Q*Yt)/G
        self.beta=beta
        rho=Yt-Ft*beta
        
        
        sigma2=0
        for i in range(len(rho)):
            sigma2=sigma2+rho[i,0]**2      
        sigma2=sigma2/self._Data.npts;
        detR=np.prod(np.diag(C)**(2.0/self._Data.npts))
        fi=sigma2*detR
        #structure with the elements
        par=Struct(beta=beta,sigma2=sigma2,C=C,Yt=Yt,Ft=Ft,G=G,R=R)
        self._par=par
        return fi

    def solve(self):
        '''
        Solves the sita minimization problem for the model's data
        
        Parameters
        ----------
        None

       
        Returns
        -------
        Optimum sita vector. 
        '''        
     
#        Minimize function using the modified powell method.
#        self._objfun: the objective function to be minimized.
#        self._sita: initial guess
#        sita_opt: updated sita that minimizes function
  
        sita_opt=fmin_powell(self._objfun,self._Data.sita,ftol=0.01,
                             maxiter=10000, disp=False)
        
##        print "optimized sita="
##        print self._sita_opt
        #Catch to be sure we always return a list of sita,
        #even if only 1 value/0-d array 
        if (len(np.atleast_1d(sita_opt)) == 1):
            sita_opt_fixed = np.array([sita_opt])
        else:
            sita_opt_fixed = sita_opt

        self._sita_opt=sita_opt_fixed
        return sita_opt_fixed
    
    def setFixed(self, sita_array):
        '''
        Set the model to use a fixed value of sita in place of solving for 
        sita numerically. Value of sita must be supplied

        Parameters
        ----------
        sita_array  np.array of sita values
        
        Returns
        -------
        Objective function value at this sita
        '''
        result = self._objfun(sita_array)                
        self._sita_opt = sita_array
        
        return result
        

    def predictor(self,xpred):
        '''
        Computes the Kriging interpolator at the provided design points(s)
        
        Parameters
        ----------
        xpred:  Numpy matrix of design points in same format of xvals for 
                constructor

       
        Returns
        -------
        Numpy matrix where column 0 is the predicted function response and 
        column 1 is the estimated mean square error.  As many rows as 
        xpred
                
        # Revised to include a preduction of Jacobian matrix of the ypred
        
        
        '''   
        f=np.matrix(np.ones(1))
        ypred=np.matrix(np.ones((xpred.shape[0],1)))
        
        ypred_d = np.matrix(np.ones((xpred.shape[0],self._Data.nvars)))   
        
        ypred_d2 = np.matrix(np.ones((xpred.shape[0],self._Data.nvars))) 
        
        ypred_hess = np.matrix(np.ones((xpred.shape[0], self._Data.nvars)))         
        
        #print ypred_hess[0].shape
        mse=np.matrix(np.ones((xpred.shape[0],1)))
        #normalize the data
        xpred=(xpred-self._mx)/self._sx
        for k in range(xpred.shape[0]):
            r=np.matrix(np.zeros((self._Data.npts,1)))
            for i in range(self._Data.npts):
                for j in range(self._Data.nvars):
                    r[i]=r[i]-self._sita_opt[j]*np.absolute(
                                    xpred[k,j]-self._Data.xpts[i,j])**2
                r[i]=np.exp(r[i])
            
            
          
        
            # compute Jacobian of r matrix
            r_j = np.matrix(np.zeros((self._Data.npts,self._Data.nvars)))
            for i in range(self._Data.npts):
                for j in range(self._Data.nvars):
                    r_j[i,j] = r[i]*(-2*self._sita_opt[j])*(xpred[k,j]-self._Data.xpts[i,j])
            
            
            r_h = np.matrix(np.zeros((self._Data.npts,self._Data.nvars)))   
            for i in range(self._Data.npts):
                for j in range(self._Data.nvars):
                    r_h[i,j] = r[i]*((-2*self._sita_opt[j])+(2*self._sita_opt[j]*(xpred[k,j]-self._Data.xpts[i,j]))**2)

            ypred[k]=self._par.beta+np.transpose(r)*np.transpose(
                  inv(self._par.C))*(self._par.Yt-self._par.Ft*self._par.beta)
                  
           
            for i in range(self._Data.nvars):                
                ypred_d[k, i] = np.transpose(r_j[:, i])*np.transpose(inv(self._par.C))*(self._par.Yt-self._par.Ft*self._par.beta)
                 
                ypred_d2[k,i] = np.transpose(r_h[:, i])*np.transpose(inv(self._par.C))*(self._par.Yt-self._par.Ft*self._par.beta)
            
            hess = np.zeros((self._Data.nvars,1,self._Data.nvars))
            
            for n in range(self._Data.nvars):
                r_hess = np.matrix(np.zeros((self._Data.npts,self._Data.nvars)))
                for i in range(self._Data.npts):
                    for j in range(self._Data.nvars):
                        r_hess[i,j] = r[i]*(-2*self._sita_opt[j])*(xpred[k,j]-self._Data.xpts[i,j])*(-2*self._sita_opt[n])*(xpred[k,n]-self._Data.xpts[i,n])
                      
                for i in range(self._Data.nvars):
                    ypred_hess[k,i] = np.transpose(r_hess[:, i])*np.transpose(inv(self._par.C))*(self._par.Yt-self._par.Ft*self._par.beta)
                hess[n] = ypred_hess[k] 
            #print hess   
            hess = np.matrix(hess) 
            

            for i in range(self._Data.nvars):
                hess[i,i] = ypred_d2[k,i]
                
            
           #unnormalize the predictor
            ypred[k]=ypred[k]*self._sy+self._my
            
            #unnormalize the predicted derivatives
            for i in range(self._Data.nvars):
                ypred_d[k, i] = ypred_d[k, i]*self._sy/self._sx[i]
                
                ypred_d2[k,i] =  ypred_d2[k,i]*self._sy/(self._sx[i])**2
            # unnormalize predicted hessian matrix
      
            for i in range(self._Data.nvars):
                for j in range(self._Data.nvars):
                    hess[i,j] = hess[i,j]*self._sy/(self._sx[i]*self._sx[j])
            #print hess    
                
            rt=inv(self._par.C)*r
            u=np.transpose(self._par.Ft)*rt-f
           
            v=1.0/self._par.G*u
            mse[k]=self._par.sigma2*self._sy**2*(1+v**2-sum(np.power(rt,2)))
            if ((mse[k] < 0) and (mse[k]/ypred > -1.0e-3)):
                #This happens when a sampled site is entered 
                mse[k] = 0
            
                
        return ypred, ypred_d, hess, mse

        
                

        
        
        
        
        
        


    

    

    
        
        
    


    
        
        

