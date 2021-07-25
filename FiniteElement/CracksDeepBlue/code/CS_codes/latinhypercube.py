
#!/usr/bin/python
# -*- coding:utf-8 -*-
#-----------------------------------------------------------------------------
# Name:        lhs.py
# Project:  Bayesian-Inference
# Purpose:     
#
# Author:      Flvio Codeo Coelho  <fccoelho@gmail.com>
#
# Created:     2008-11-26
# Copyright:   (c) 2008 by the Author
# Licence:     GPL
#-----------------------------------------------------------------------------
__docformat__ = "restructuredtext en"
#from pylab import plot, figure,hist,show, savefig, legend
import scipy.stats as stats
import numpy

def lhsFromSample(sample,siz=100):
    """
    Latin Hypercube Sample from a set of values

    :Parameters:
        - `sample`: list, tuple of array
        - `siz`: Number or shape tuple for the output sample
    """
    if not isinstance(sample, (list,tuple,numpy.ndarray)):
        raise TypeError('sample is not a list, tuple or numpy vector')
    n = siz
    if isinstance(siz,(tuple,list)):
        n=numpy.product(siz)
    perc = numpy.arange(0,100.,100./n)
    numpy.random.shuffle(perc)
    smp = [stats.uniform(i,100./n).rvs() for i in perc]
    v = numpy.array([stats.scoreatpercentile(sample,p) for p in smp])
    if isinstance(siz,(tuple,list)):
        v.shape = siz
    return v

def lhsFromDensity(kde,siz=100):
    '''
    LHS sampling from a variables Kernel density estimate.

    :Parameters:
        - `kde`: scipy.stats.kde.gaussian_kde object
        - `siz`: Number or shape tuple for the output sample
    '''
    if not isinstance(kde,scipy.stats.kde.gaussian_kde):
        raise TypeError("kde is not a density object")
    if isinstance(siz,(tuple,list)):
        n=numpy.product(siz)
    s = kde.resample(n)
    v = lhsFromSample(s,n)
    if isinstance(siz,(tuple,list)):
        v.shape = siz
    return v


def lhs(dist, parms, siz=100):
    '''
    Latin Hypercube sampling of any distrbution.
    dist is is a scipy.stats random number generator 
    such as stats.norm, stats.beta, etc
    parms is a tuple with the parameters needed for 
    the specified distribution.

    :Parameters:
        - `dist`: random number generator from scipy.stats module.
        - `parms`: tuple of parameters as required for dist.
        - `siz` :number or shape tuple for the output sample
    '''
    if not isinstance(dist, (stats.rv_discrete,stats.rv_continuous)):
        raise TypeError('dist is not a scipy.stats distribution object')
    n=siz
    if isinstance(siz,(tuple,list)):
        n=numpy.product(siz)
    #force type to float for sage compatibility
    parms = tuple([float(k) for k in parms])
    perc = numpy.arange(0.,1.,1./n)
    numpy.random.shuffle(perc)
    smp = [stats.uniform(i,1./n).rvs() for i in perc]
    v = dist(*parms).ppf(smp)
    
    if isinstance(siz,(tuple,list)):
        v.shape = siz
    return v
            
if __name__=='__main__':
    dist = stats.uniform
    #dist = stats.beta
    pars = (0,1)
    #pars = (1,5) #beta
    c=lhs(dist, pars,(20,5))
    print c
    #hist(c,normed=1, label='LHS sample')
    n = dist(*pars).rvs(size=20)
    #hist(n.ravel(),facecolor='r',alpha =0.3,normed=1, label='Regular sample')
    #plot(numpy.arange(min(min(c),min(n)),max(max(c),max(n)),.1),dist(*pars).pdf(numpy.arange(min(min(c),min(n)),max(max(c),max(n)),.1)),label='PDF')
    #legend()
    #savefig('lhs.png',dpi=400)
    #show()
    

#TODO: Add correlated multiple lhs sampling