#*****************************************************************
#Copyright (c) 2016 Regents of the University of Michigan
#Author: Yan Liu and Matt Collette
#Contact: mdcoll@umich.edu
#See LICENSE files for details of license and code 
#*****************************************************************
from __future__ import division
import numpy as np
import copy
from random import Random
import sys
import MSDL_SOGA as soga

__author__ = 'Dylan Temple'


class Selection(object):
    '''Various selection operators for evolutionary operators'''
    def __init__(self, constraint=None,rng=None,operator='two-pass-tourney'):
        '''
        Initialize the random number generator for the selection operator
        and the constraint instance to be used in the comparisons
        
        Parameters
        ----------
        constraint: An instance of constraint_handler class for use in comparison
                        of individuals.  If None is given no constraint data will be used
                        in the comparison (default=None).
        rng:        An instance of a psuedo-random number generator to be used in the
                        stochastic aspects of the operator.  If None is given then a new
                        instance of Random() will be generated to use (default=None)
        operator:   The type of selection operator for the 'evaluate' method to use if this
                        this class is being called via its initial instanciation. (default='two-pass-tourney')          
        '''
        if rng:
            self.rng = rng
        else:
            self.rng = Random()
        
        self.constraint = constraint
        self.operator = operator
    
    def _compare(self, ind1, ind2):
        '''
        Compare two individuals using the prescribed constraint handling
        technique.
        
        Parameters
        ----------
        ind1:   The first individual to be compared
        ind2:   The second individual to be compared
        
        Returns
        --------
        winner: The better of the two individuals
        '''
        #Get scaled objective information
        obj_1 = ind1.scaled_objective()
        obj_2 = ind2.scaled_objective()
        
        #Determine better solution
        if obj_1 < obj_2:
            winner = ind1
        else:
            winner = ind2
    
        return winner
            
    def _two_pass_tournament(self, population):
        '''
        A two pass torunament type selection
        
        Parameters
        ----------
        population: The population of individuals to select from
        
        Returns
        --------
        survivors:  The selected population
        '''
        survivors = []
        
        #Perform the tournament selection twice.  Shuffle
        #the population each time
        for i in range(0,2):
            pop_copy = copy.copy(population)
                
            self.rng.shuffle(pop_copy)
            for j in range(0, (len(pop_copy)//2)*self.pop_percent):
                #Select 2 individuals and compare
                ind1 = pop_copy.pop()
                ind2 = pop_copy.pop()
                winner = self._compare(ind1, ind2)
                survivors.append(winner)
        return survivors
    
    def _roulette_wheel(self, population):
        '''
        A roulette wheel type selection
        
        Parameters
        ----------
        population: The population of individuals to select from
        
        Returns
        --------
        survivors:  The selected population
        '''
        survivors = []
        objectives = []
        #Get scaled objective vector for entire population
        for ind in population:
            objectives.append(ind.scaled_objective())
        
        #Scale objectives to be used in roulette wheel selection (simulate a 
        #maximization problem)
        scaled_objectives = []
        max_obj = np.max(objectives)
        for obj in objectives:
            scaled_objectives.append(max_obj - obj)
        
        #Create roulette bins for selection
        obj_sum = sum(scaled_objectives)
        obj_contribution = []
        for obj in scaled_objectives:
            obj_contribution.append(obj/obj_sum)
        bins = np.cumsum(obj_contribution)
        #Spin wheel N times and populate survivor list
        for i in range(len(population)*self.pop_percent):
            spin = self.rng.random()
            index = next(n for n, j in enumerate(bins) if j > spin)
            survivors.append(population[index])

        return survivors
        
    
    def evaluate(self, population, pop_percent=1):
        '''
        Selects a segment of a population based on a method described in the 
        initialization of the class.
        
        Parameters
        ----------
        population:     A list of individual classes making up the population to be
                            selected from.
        pop_percent:    The percent size of the new population compared to the old (default=1.0)
        
        Returns
        -------
        new_pop:    The selected portion of the population chosen via the strategy
                        chosen when the selection class was initialized.
        '''
        self.pop_percent = pop_percent
        
        #Select a population
        if self.operator.lower() == 'two-pass-tourney':
            new_pop = self._two_pass_tournament(population)
        elif self.operator.lower() == 'roulette-wheel':
            new_pop = self._roulette_wheel(population)
        else:
            sys.exit("\nThe selection class takes operators of 'two-pass-tourney' or 'roulette-wheel'.  Got '"+str(self.operator)+"' instead.")
        return new_pop
            
class Constraint_Handler(object):
    '''Various constraint handling operators'''
    
    def __init__(self, operator='quadratic-penalty'):
        '''
        Initialize the type of constraint handling technique that will be associated
        with this instance.
        
        Parameters
        ----------
        operator:   The type of constraint handling approach to associate with this
                        class instance. (default='quadratic-penalty)
        '''
        self.operator = operator
    
    def _quadratic_penalty(self, ind):
        '''
        A quadratic exterior penalty method constraint handler
        
        Parameters
        -----------
        ind:    The individual to return the constraint value for
        
        Returns
        --------
        constr: The scaled constraint violation for the individual
        '''
        
        #Get the sum of the squares of the constraint violation vector
        penalty = 0.
        violation = ind.violation()
        for i in range(len(violation)):
            penalty += violation[i]**2
        
        #Apply the scaling parameter
        constr = (ind.constr_param/2) * penalty
        return constr
        
    def _log_barrier(self, ind):
        '''
        A logarithmic interior penalty method constraint handler
        
        Parameters
        -----------
        ind:    The individual to return the constraint value for
        
        Returns
        --------
        constr: The scaled constraint violation for the individual
        '''
        #Get the sum of the log of the constraint violation vector
        penalty = 0.
        violation = ind.violation()
        for i in range(len(violation)):
            penalty += np.log(violation[i])
        
        #Apply the negative of the scaling parameter
        constr = -1 * ind.constr_param * penalty
        return constr
    
    def _inverse_barrier(self, ind):
        '''
        An inverse interior penalty method constraint handler
        
        Parameters
        -----------
        ind:    The individual to return the constraint value for
        rho:    The scaling parameter for the penalty function
        
        Returns
        --------
        constr: The scaled constraint violation for the individual
        '''     
        #Get the sum of the inverse of the constraint violation vector
        penalty = 0.
        violation = ind.violation()
        for i in range(len(violation)):
            if violation[i] > 0.:
                penalty += violation[i]**-1
        
        #Apply the scaling parameter
        constr = ind.constr_param * penalty
        return constr
    
    def _feasible_penalty(self, ind):
        '''
        A feasible domain penalty function proposed in Deb 2001.  
        
        Parameters
        ----------
        ind:    The individual to return the scaled value for
        fmax:   The maximum objective function value in the population
        
        Returns
        --------
        constr: The scaled constraint violation for the individual
        '''
        #Get constraint violation and sum
        violation_sum = sum(ind.violation())
        
        if violation_sum > 0.:
            constr = violation_sum + ind.constr_param
        else:
            constr = 0.
        
        return constr
        
    def evaluate(self, ind):
        '''
        Evaluate the constraint value based on the constraint handling technique
        defined in the initialization of the class.
        
        Parameters
        ----------
        ind:    Individual to get the constraint value for
        
        Returns
        -------
        constr: The constraint penalty term for the individual
        '''
        if self.operator.lower() == 'quadratic-penalty':
            constr = self._quadratic_penalty(ind)
        elif self.operator.lower() == 'log-barrier':
            constr = self._log_barrier(ind)
        elif self.operator.lower() == 'inverse-barrier':
            constr = self._inverse_barrier(ind)
        elif self.operator.lower() == 'feasible-penalty':
            constr = self._feasible_penalty(ind)
        else:
            sys.exit("A constraint_handler class must be initialized with an operator of 'external-penalty', \
                        'log-barrier', 'inverse barrier', or 'feasible-penalty'.  Got '"+str(self.operator)+"' instead")
        
        return constr


class Crossover(object):
    '''Various crossover operators for evolutionary algorithms'''
    def __init__(self, rng=None, operator='simulated-binary', cpercent=0.0):
        '''
        Initialize the crossover classs
        
        Parameters
        ---------
        rng:        An instance of a pseudo-random number generator.  If None is supplied
                        then a new instance of Random() will be used. (default=None)
        operator:   The type of crossover operator approach to associate with this
                        class instance. (default='quadratic-penalty)
        cpercent:   The percent chance for crossover to occur.  (default=0.9)
        '''
        if rng:
            self.rng = rng
        else:
            self.rng = Random()
        
        self.operator=operator
        self.cpercent = cpercent
    
    def _simulated_binary(self, parent1, parent2, eta=4.0):
        '''
        Creates two child genes from two parent genes via simulated binary
        crossover.  This is for real coded continuous genes only.
        
        Parameters
        ---------
        parent1:    The first parent gene
        parent2:    The second parent gene
        eta:        The crossover exponent (default=4.0)
        
        Returns
        --------
        child1:     The first child gene
        child2:     The second child gene
        '''
        pc = self.rng.random()
        if pc >= self.cpercent:
            #If crossover doesn't happen return parent genes
            child1 = parent1
            child2 = parent2
        else:
            #Perform crossover
            u = self.rng.random()
            if u < 0.5:
                Bq = (2*u)**(1/(eta+1))
            else:
                Bq = (1/(2*(1-u)))**(1/(1+eta)) 
            child1=0.5*((1.0+Bq)*parent1.value+(1.0-Bq)*parent2.value)
            child2=0.5*((1.0-Bq)*parent1.value+(1.0+Bq)*parent2.value)   
            if child1 > parent1.up_bound:
                child1 = parent1.up_bound
            elif child1 < parent1.low_bound:
                child1 = parent1.low_bound
            if child2 > parent2.up_bound:
                child2 = parent2.up_bound
            elif child2 < parent2.low_bound:
                child2 = parent2.low_bound
            child1 = soga.Gene(child1, parent1.up_bound, parent1.low_bound, parent1.gtype, parent1.encoding)
            child2 = soga.Gene(child2, parent2.up_bound, parent2.low_bound, parent2.gtype, parent2.encoding)
        return child1, child2
    
    def _single_point(self, parent1, parent2):
        '''
        Creates two child child genes from two parent genes using single
        point crossover.  This is for binary coded genes only.
        
        Parameters
        ----------
        parent1:    The first parent gene
        parent2:    The second parent gene
        
        Returns
        --------
        child1:     The first child gene
        child2:     The second child gene
        '''
        pc = self.rng.random()
        if pc >= self.cpercent:
            #If crossover doesn't happen return parent genes
            child1 = parent1
            child2 = parent2
        else:       
            #Perform crossover
            cross_point = self.rng.randint(1,len(parent1.value)-1)
            child1 = parent1.value[:cross_point] + parent2.value[cross_point:]
            child2 = parent2.value[:cross_point] + parent1.value[cross_point:]
            child1 = soga.Gene(child1, parent1.up_bound, parent1.low_bound, parent1.gtype, parent1.encoding)
            child2 = soga.Gene(child2, parent2.up_bound, parent2.low_bound, parent2.gtype, parent2.encoding)
            
        
        return child1, child2
    
    def _double_point(self, parent1, parent2):
        '''
        Creates two child child genes from two parent genes using double
        point crossover.  This is for binary coded genes only.
        
        Parameters
        ----------
        parent1:    The first parent gene
        parent2:    The second parent gene
        
        Returns
        --------
        child1:     The first child gene
        child2:     The second child gene
        '''
        pc = self.rng.random()
        if pc >= self.cpercent:
            #If crossover doesn't happen return parent genes
            child1 = parent1
            child2 = parent2
        else:       
            #Perform crossover
            cross_points = (self.rng.randint(1,len(parent2.value)-1), self.rng.randint(1,len(parent2.value)-1))
            cross_points = sorted(cross_points)
            child1 = parent1.value[:cross_points[0]] + parent2.value[cross_points[0]:cross_points[1]] + parent1.value[cross_points[1]:]
            child2 = parent2.value[:cross_points[0]] + parent1.value[cross_points[0]:cross_points[1]] + parent2.value[cross_points[1]:]
            child1 = soga.Gene(child1, parent1.up_bound, parent1.low_bound, parent1.gtype, parent1.encoding)
            child2 = soga.Gene(child2, parent2.up_bound, parent2.low_bound, parent2.gtype, parent2.encoding)

        return child1, child2
    
    def _multi_point(self, parent1, parent2, points):
        '''
        Creates two child child genes from two parent genes using multi
        point crossover.  This is for binary coded genes only.
        
        Parameters
        ----------
        parent1:    The first parent gene
        parent2:    The second parent gene
        multi:      The number of cross over points
        
        Returns
        --------
        child1:     The first child gene
        child2:     The second child gene
        '''        
        
        if points == 1 or points == 2:
            print 'Multi point cross over operator being used for single or double point crossover is inefficient.  Use, \
                    single-point or double-point crossover for better performance'
        if type(points).__name__ != 'int':
            syx.exit("Multi Point crossover needs an integer value for the number of points.  Got '"+str(points)+"' instead.")
                    
        pc = self.rng.random()
        if pc >= self.cpercent:
            #If crossover doesn't happen return parent genes
            child1 = parent1
            child2 = parent2
        else:       
            #Perform crossover
            cross_points = []
            for i in range(points):
                cross_points.append(self.rng.randint(1,len(parent1.value)-1))
            cross_points = sorted(cross_points)
            parents = (parent1, parent2)
            children = []
            parent = False
            for i in range(0,2):
                child = []
                for j in range(points+1):
                    if j == 0:
                        child += parents[parent].value[:cross_points[j]]
                        parent = not parent
                    elif j == points:
                        child += parents[parent].value[cross_points[j-1]:]
                    else:
                        child += parents[parent].value[cross_points[j-1]:cross_points[j]]
                        parent = not parent
                children.append(child)
                parent = True
            
            child1, child2 = children[0], children[1]
            child1 = soga.Gene(child1, parent1.up_bound, parent1.low_bound, parent1.gtype, parent1.encoding)
            child2 = soga.Gene(child2, parent2.up_bound, parent2.low_bound, parent2.gtype, parent2.encoding)

        return child1, child2
    
    def evaluate(self, parent1, parent2, param=None):
        '''
        Performs the cross over operation selected upon initialization of the class
        for two genes.
        
        Parameters
        ----------
        parent1:    The first parent gene
        parent2:    The second parent gene
        param:      A parameter for the cross over operator, i.e. eta for sbx
                    or the number of points for multi point.(default=None)
        
        Returns
        --------
        child1:     The first child gene
        child2:     The second child gene
        ''' 
        
        if self.operator.lower() == 'simulated-binary':
            child1, child2 = self._simulated_binary(parent1, parent2, param)
        elif self.operator.lower() == 'single-point':
            child1, child2 = self._single_point(parent1, parent2)
        elif self.operator.lower() == 'double-point':
            child1, child2 = self._double_point(parent1, parent2)
        elif self.operator.lower() == 'multi-point':
            child1, child2 = self._multi_point(parent1, parent2, param)
        else:
            sys.exit("The availavle crossover operators are 'simulated-binary', 'single-point', \
                        'double-point', or 'multi-point'. Got '"+str(self.operator)+"' instead.")
        
        return child1, child2
        
class Mutation(object):
    '''Various mutation operators for evolutionary optimizers'''
    def __init__(self, rng=None, operator='nsgaII', mpercent=0.005):
        '''
        Initialize the mutation operator class
        
        Parameters
        ---------
        rng:        An instance of a pseudo-random number generator.  If None is supplied
                        then a new instance of Random() will be used. (default=None)
        operator:   The type of mutation operator approach to associate with this
                        class instance. (default='quadratic-penalty)
        mpercent:   The percent chance that mutation occurs
        '''
        if rng:
            self.rng = rng
        else:
            self.rng = Random()
        
        self.operator = operator
        self.mpercent = mpercent
    
    def _nsgaII(self, gene, mut_x=2.0):
        '''
        Mutate a gene using the method in Deb's NSGA-II algorithm (2002).  This
        operator is for real coded continuous genes.
        
        -----------
        gene:   The gene to be mutated
        
        Returns
        --------
        gene:       The mutated gene
        mutFlag:    Returns True if mutation happens, False otherwise
        '''
        mutFlag = False
        max_mutation = gene.up_bound - gene.low_bound
        pm = self.rng.random()
        u = self.rng.random()
        if pm >= 1.00-self.mpercent:
            mutFlag = True
            if u < 0.5:
                delta = (2*u)**(1/(1+mut_x))-1
            else:
                delta = 1-(2*(1-u))**(1/(mut_x+1))
        else:
            delta = 0      
            
        gene.value += delta*max_mutation
        if gene.value > gene.up_bound:
            gene.value = gene.up_bound
        elif gene.value < gene.low_bound:
            gene.value = gene.low_bound
        return gene
    
    def _binary(self, gene):
        '''
        Mutate a binary gene by simple digit flipping.  This operator is only
        for use with binary coded variables.
        
        Parameters
        -----------
        gene:   The gene to be mutated
        
        Returns
        --------
        gene:   The mutated gene
        '''
        for i in range(len(gene.value)):
            pm = self.rng.random()
            if pm >= self.mpercent:
                gene.value[i] = gene.value[i]
            else:
                gene.value[i] = int(not gene.value[i])
        
        return gene
    
    def evaluate(self, gene, param=None):
        '''
        Mutates a gene using the prescribed mutation operator.
        
        Parameters
        -----------
        gene:   The gene to be mutated
        param:  The parameter for mutation, i.e. the mutation exponent in
                    the NSGAII operator.  (default=None)
        
        Returns
        --------
        gene:   The mutated gene
        '''
        
        if self.operator.lower() == 'nsgaii':
            gene = self._nsgaII(gene, param)
        elif self.operator.lower() == 'binary':
            gene = self._binary(gene)
        else:
            sys.exit("The available operators for mutation are 'nsgaII' and 'binary'.  Got, \
                        '"+str(self.operator)+"' instead")
        
        return gene
        
        
            
        
        
        
                    
            
        
                    
                
        
    
        
        
    
    
    
    

                        
        
        
        
        
            
    
    
    
    