#*****************************************************************
#Copyright (c) 2016 Regents of the University of Michigan
#Author: Yan Liu and Matt Collette
#Contact: mdcoll@umich.edu
#See LICENSE files for details of license and code 
#*****************************************************************

import copy
import numpy as np
import math

class Comparison:
    '''
    Class for comparing different pareto fronts obtain by implemting different 
    optimizing methods.
    
    Size, Coverage, and Non-Uniformity taken from  from:
    Comparison of Multi-Objective Genetic Algorithms in Optimizing Q-Law 
    Low-Thrust Orbit Transfers written by Seungon Lee, Paul v on Almen, 
    Wolfgang Fink, Anastassios E. Petropoulos, and Richard J. Terrile
    
    Other methods are developed in-house
    '''
    '''
    def __init__ (self):
    '''    
    def Size (self,obj1_1,obj2_1,obj3_1,obj1_2,obj2_2,obj3_2):
        '''
        Size of the dominated space, measure how much a non-dominated set span.
        The larger the value the better the solution
        '''
        max1_1=max(obj1_1)
        max2_1=max(obj2_1)
        max3_1=max(obj3_1)
        print max1_1
        print max2_1
        print max3_1
        max1_2=max(obj1_2)
        max2_2=max(obj2_2)
        max3_2=max(obj3_2)
        print max1_2
        print max2_2
        print max3_2
        if max1_1 >= max1_2 and max2_1 >= max2_2 and max3_1 >= max3_2:
            value_1 = 1
            value_2 = -1
        elif max1_1 < max1_2 and max2_1 < max2_2 and max3_1 < max3_2:
            value_1 = -1
            value_2 = 1
        else:
            value_1 = value_2 =0
        return value_1,value_2
    def Coverage (self,obj1_1,obj2_1,obj3_1,obj1_2,obj2_2,obj3_2,tolerance):
        '''
        Coverage measures two Pareto fronts obtained from different methods. 
        The higher values means the Pareto front dominate the other Pareto front
        better.
        '''
        method1=[]
        method2=[]
        count1 = 0
        count2 = 0
        for i in range(len(obj1_1)):
            obj=[obj1_1[i],obj2_1[i],obj3_1[i]]
            method1.append(obj)
        for i in range(len(obj1_2)):
            obj=[obj1_2[i],obj2_2[i],obj3_2[i]]
            method2.append(obj)
        copymethod2 = copy.copy(method2)
        for k in range(len(obj1_1)):
            #print "number:", k
            o1=method1[k]
            for o2 in copymethod2:
                flag = 1
                for i in range(0,len(o1)):
                    if o1[i]>o2[i]:
                        flag=0
                    elif o1[i]<o2[i]:
                        if abs((o1[i]-o2[i])/o2[i])>=tolerance:
                            flag*=2
                        else:
                            flag=0
                if flag>1:
                    #method1 dominate method2
                    copymethod2.remove(o2)
                    count1 +=1
                    #print "dom"
        print "the number of method 1 dominate method2", count1
        copymethod1 = copy.copy(method1)
        for k in range(len(obj1_2)):
            #print "number:", k
            o2=method2[k]
            for o1 in copymethod1:
                flag = 1
                for i in range(0,len(o1)):
                    if o2[i]>o1[i]:
                        flag=0
                    elif o2[i]<o1[i]:
                        if abs((o2[i]-o1[i])/o1[i])>=tolerance:
                            flag*=2
                        else:
                            flag=0
                if flag>1:
                    #method2 dominate method1
                    copymethod1.remove(o1)
                    count2 +=1
                    #print "dom"
                    #print "o1",o1
                    #print "o2",o2
        print "the number of method 2 dominate method1", count2
        return float(count1)/len(method2),float(count2)/len(method1)
            
    def NonUniformity (self,crowd1,crowd2):
        '''
        Non-uniformity meassure the distribution of the final Pareto fronts 
        obtained by different methods. The lower of the value means the Pareto
        front is more evenly distributed and is more preferable.
        '''
        count1 = 0
        sum1=0
        for i in crowd1:
            if i >= 10000:
                crowd1.remove(i)
                count1 +=1
            else:
                sum1 +=i
        for j in range (0,2):
            for i in crowd1:
                if i >= 10000:
                    crowd1.remove(i)
        #print "# of boundaries", count1
        #print sum(crowd1)
        average1 = sum1 / len(crowd1)
        #print "average1:", average1
        #print "average1:", average1
        deviation1=[]
        for i in crowd1:
            dev = ((i / average1)-1)**2
            deviation1.append(dev)
        dist1 = float(sum(deviation1))/ len(crowd1)
        sum2=0
        count2 = 0
        for i in crowd2:
            if i >= 10000:
                crowd2.remove(i)
                count2 +=1
            else:
                sum2 +=i
        for i in crowd2:
            if i>=10000:
                crowd2.remove(i)
        #print "# of boundaries", count2
        average2 = sum2 / len(crowd2)
        #print "average2",average2
        #print "average2:", average2
        deviation2=[]
        for i in crowd2:
            dev = ((i / average2)-1)**2
            deviation2.append(dev)
        dist2 = float(sum(deviation2))/ len(crowd2)
        return dist1,dist2

    def NormalizedDistance (self,refPop,assessPop):
        '''
        Calculates a metric of how close the members of assessPop lie
        to refPop
                 
        Developed by Matthew Collette.  refPop and assessPop are assumed
        to be lists of objective function tuples/lists. For each member in 
        assessPop the normalized distance to all members in refPop is calculated
        and then the lowest distance is kept for each member.  The mean 
        and standard deviation of these normalized distances is then returned.
        
        Parameters
        ----------
        refPop:  List-like, each row represents the objective functions of 
                 a member of the population. Reference population that distances
                 are measured too.
                 
        assessPop: List-like, each row represents the objective functions of
                  a member of the population. This is the population to be 
                  assessed
                  
        Returns
        -------
        Tuple of the mean minimum distances and the standard deviation of 
        this distance (sample n-1 based st. dev)

       '''
        #Logical result if no assess pop - if Pareto front not existing for
        #some reason
        if len(assessPop) == 0:
            return (0,0)
        #Distances from a member in assessPop to all members of refPop
        distances = np.empty((len(refPop)))    
        #minimum distanes from each member in assessPop to any member of refPop
        min_distances = np.empty((len(assessPop)))
        
        #Tolerance for zero values
        TOL = 1.e-14
         
        #Number of objective functions
        num_obj = len(refPop[0])
        i = 0
        for ap in assessPop:
            j = 0
            for rp in refPop:
                total = 0
                for k in xrange(0, num_obj):
                    if math.fabs(rp[k]) > TOL:
                        total = total + ((rp[k] - ap[k])/rp[k])**2.0
                    else:
                        total = total + ap[k]**2.0
                distances[j] = (total)**0.5
                j = j + 1
            min_distances[i] = np.amin(distances)
            i = i + 1
        return (np.mean(min_distances), np.std(min_distances, ddof=1))

         
         
        
         