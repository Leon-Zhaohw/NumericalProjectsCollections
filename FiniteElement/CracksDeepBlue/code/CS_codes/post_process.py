#*****************************************************************
#Copyright (c) 2016 Regents of the University of Michigan
#Author: Yan Liu and Matt Collette
#Contact: mdcoll@umich.edu
#See LICENSE files for details of license and code 
#*****************************************************************

import subprocess
import sqlite3
import numpy as np
import matplotlib.pyplot as plt
import math
#import methods.regression.kriging as kriging
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm



def isDom(point_1, point_2):
    '''Determines if one point dominates another for any number of obj. functs
    
    Parameters
    ----------
    point_1:    List-like of objective function values 1...n for point 1

    point_2:    List-like of objective function values 1...n for point 2

    Returns
    -------
    -1 if point 1 dominates point 2, 1 if 2 dominates 1, and 0 if points are 
    non-dominated    
    '''
    oneDomTrue = 1
    twoDomTrue = 1
    for (pt1, pt2) in zip(point_1, point_2):
        if pt1 < pt2:
            oneDomTrue *= 2
            twoDomTrue = 0
        if pt2 < pt1:
            oneDomTrue = 0
            twoDomTrue *=2
    
    assert (not(oneDomTrue==2 and twoDomTrue==2)), "failed domination status"
    
    #Figure out if they are dominated or not
    if (oneDomTrue==twoDomTrue):
        return  0
    if (oneDomTrue > 1):
        return -1
    if (twoDomTrue > 1):
        return 1
    
    assert (false), "failed domination status"
    return


def combinePareto(input_list):
    '''Pareto-only sort on one or more pt list
    
    Parameters
    ----------
    input_list:  triply-nested list, where each first-level sublist is in the 
                 format returned by getFront on NSGA post-processor 
    
    Returns
    --------
    Tuple of two lists, one is a list of list of objective function values,
    each list corrosponding to an objective function, not an individual,
    suitable for plotting of the non-dominated points in the input, and 
    the second is a list of two-item tuples listing input set data number
    (starting with 0) and individual ID number of the corrosponding points
    '''
    
    #This is just an ugly O(n^2) loop could be much faster if list built 
    pareto_list = []
    index = 0
    for dataset in input_list:
        for pt in dataset:
            okToAdd = True
            for compareset in input_list:
                for pt_compare in compareset:
                    if (isDom(pt[1:], pt_compare[1:])==1):
                        okToAdd = False            
            if (okToAdd):
                pareto_list.append((pt, index))
        index += 1
    
    #Now go through the output list and make the correct formating
    num_obj = len(input_list[0][0]) - 1
    output_objs = []
    for k in range(num_obj):
        output_objs.append([])

    data_set_and_ids = []
    for pt in pareto_list:
        for count in range(1, num_obj +1):
            output_objs[count-1].append(pt[0][count])
        data_set_and_ids.append((pt[1], pt[0][0]))
    
    return (output_objs, data_set_and_ids)

                
class NSGA_PostProcess:
    '''
    Class for post-processing NSGA-II data stored in SQL databases. 
    
    '''
    def __init__(self, database):
        '''
        Initiates the object, setting up a connection to the provided 
        database
        '''
        self.conn=sqlite3.connect(database)
        self.cur=self.conn.cursor()
    
    
    def addFrontTo2DPlot(self, front_list, ax= None, index_x = 1, index_y = 2,
                       format_string = 'ro', forced_axis=None):
        '''Utility function to add a front to a 2-D plot from the format 
        returned by getFront
        
        Parameters
        ----------
        front_list:     A list of points in the front in the format returned by
                        getFront
                        
        ax:             Current axis, if none uses plt.gca() to find
        
        index_x:        List index of the front to be plotted on the x-axis.
                        Note that this is not the NSGA-II objective number
                        rather the list index where the front is returned 
                        front getFront.  Defaults to 1, the first front after 
                        ID number in the list
        
        index_y:        List of index of y-axis front, as above, defaults to 2
        
        format_string:  matplotlib-format string for points, defaults to 'ro'
       
        Returns
        -------
        None
        '''
        #First, reformat the data from lists-of-lists to np arrays.
        num_pts = len(front_list)
        x_vals = np.empty((num_pts,1))
        y_vals = np.empty((num_pts,1))
        kount = 0
        for pt in front_list:
            x_vals[kount] = pt[index_x]
            y_vals[kount] = pt[index_y]
            kount += 1
        
        #The do the plot 
        if ax is None:
            ax = plt.gca()
        ax.plot( x_vals, y_vals, format_string)
        if forced_axis:
            plt.axis(forced_axis)
        
        return 
    
    
    
    def getFront(self, gen_number, front_number, obj_functions, VFO=None):
        '''Returns the individual ids and selected objective function 
        values for a given front and optionally VFO status
        
        Gets the objective function values for a given front at a given
        generation number from the optimizer. Not all objective function
        values must be returned, a list of objective function values can
        be provided, and only those objective function values will be 
        returned.  generation number and front number must be scalar.
        Optionally, a tuple may be specifie for VFO.  The tuple contains
        and objective function number and a true/false e.g. (2, False)
        If specified, function will only return members of the front 
        where the objective function has (True) or has not (False) been
        updated by the VFO model.  A tolerance of 1e-10 is used on the VFO
        scale to differentiate betwen 1.0 (no scale) and scaling values
        
        Parameters
        ----------
        gen_number: scalar generation number to pull the data from
        
        front_number: scalar front number within that generatino to pull data 
                      from. Remember pareto front has index 0!
    
        obj_functions: List-like listing of objective function IDs to pull
                       corrosponding objective function for individuals
                       
        VFO: Optional, tuple of objective function ID and True/False.  If
             provided, method will only return points in front where VFO
             has (true) or has not (false) updated the objective function
        
        Returns
        -------
        List-of-lists where the first entry in each sub-list is the individual
        ID, followed by the objective function values corrosponding to the
        passed obj_functions IDs
        
        '''
        VFOTOL = 1.e-10
        #Get all individuals in the requested front at the reguested 
        #generatrion        
        #Table create command for reference: 
        #'''CREATE TABLE generation (GenNum,IndID,Front
        select_tuple = (gen_number, front_number)
        self.cur.execute('select * from generation where GenNum = ? AND ' + \
                            'FRONT = ?', select_tuple)
        #Build up a list of individuals IDS from the previous SQL statement
        ID_list = []
        for row in self.cur:
            ID_list.append(row[1]) #ID is second column in table
        
        #Go through a build objective function list 
        return_list =[]
        for indiv in ID_list:
            #Flag only used if VFO filter is active
            includePt = True
            curr_list = [indiv]
            #Go through all of the requested objective functions and retrieve
            #values
            #Table function for reference 
            #'''CREATE TABLE objfun (IndID,ObjID,Generation,Value)''' 
            for obj in obj_functions:
                sel_tuple = (indiv, gen_number, obj)
                self.cur.execute('select * from objfun where IndID= ? AND ' + \
                'Generation = ? AND ObjID = ?',
                                 sel_tuple)
                for row in self.cur:        #Should be only one
                    curr_list.append(row[3])
                    #Now check if VFO filter is on and applies to this obj
                    if (VFO and (obj == VFO[0])):
                        #Check the two cases that lead to rejection
                        #VFO=True but scale == 1.0 (no scale) or
                        #VFO=False and scale > 
                        if (VFO[1] and math.fabs(row[4]-1.) < VFOTOL):
                            includePt = False
                        elif ((VFO[1]==False) and (math.fabs(row[4] -1.) 
                                                            > VFOTOL)):
                            includePt = False
            #If no filter, or passed filter, then add to return list
            if includePt:                            
                return_list.append(curr_list)
            
        return return_list
    
    
    def getIndVariables(self, pts_list):
        '''
        Finds and returns a list-of-list of the chromosome values corrosponding
        to a given Pareto front
        
        Parameters
        ----------
        pts_list:   list-like of list of points, each is assumed to start with
                    an ID value, then optionally contain objective functions
                    In format returned by getFront
        
        Returns
        -------
        List-of-lists, where each sub-list is the chromosome of the
        corresponding individual in the passed pts_list
        '''
        #Be sure to handle empty call correctly
        if len(pts_list) == 0:
            return []
        
        #Unfortunately, there is no listing of the number of genes explicitly
        #So figure out how long the chromosome is first
        sel_tuple = (pts_list[0][0],)
        self.cur.execute('select * from chromosomes where IndID=?', sel_tuple)
        count = 0
        #Sqlite uses sequential retrivial, so need to manually count rows!
        for row in self.cur:
            count += 1
        num_genes = count        
    
    
    
        #Go through and get all the genes
        return_list = []
        
        for indiv in pts_list:
            chromosome = []
            #Go through all gene values one at a time -genes are 1 indexed
            for gene in range(1, num_genes+1):
                sel_tuple = (indiv[0], gene)
                self.cur.execute(
                 'select * from chromosomes where IndID=? and ChromID=?', 
                  sel_tuple)
                #This should only fire once - as IndID/ChromID is unique
                for row in self.cur:
                    gene = row[2]
                #Add to this individual chromosome
                chromosome.append(gene)
            #Add individual to return list
            return_list.append(chromosome)
            
        return return_list
    
    def reduceFront(self, pts_list):
        '''
        Strips leading ID term from a front list that has been returned by
        getFront
        
        Parameters
        ----------
        pts_list: list-like of list of points, each point is  assumed to have an 
                  ID value, followed by a number of objective function values
                  in format returneded by getFront. 
        
        Returns
        -------
        List-of-lists where each entry is objective function values.  Points
        are in the same order, 
        
        '''
        ret_list = []
        if len(pts_list) > 0:
            size = len(pts_list[0])
            for pt in pts_list:
                ret_list.append(pt[1:size])
        
        return ret_list
    
    def frontPlotFormat(self,gen_number,front_number, obj_functions,VFO=None):    
        '''Returns a list-of-list by objective function number for plotting
        
        Gets the objective function values for a given front at a given
        generation number from the optimizer. Not all objective function
        values must be returned, a list of objective function values can
        be provided, and only those objective function values will be 
        returned.  generation number and front number must be scalar.
        Optionally, a tuple may be specifie for VFO.  The tuple contains
        and objective function number and a true/false e.g. (2, False)
        If specified, function will only return members of the front 
        where the objective function has (True) or has not (False) been
        updated by the VFO model.  A tolerance of 1e-10 is used on the VFO
        scale to differentiate betwen 1.0 (no scale) and scaling values
        
        Parameters
        ----------
        gen_number: scalar generation number to pull the data from
        
        front_number: scalar front number within that generatino to pull data 
                      from. Remember pareto front has index 0!
    
        obj_functions: List-like listing of objective function IDs to pull
                       corrosponding objective function for individuals
                       
        VFO: Optional, tuple of objective function ID and True/False.  If
             provided, method will only return points in front where VFO
             has (true) or has not (false) updated the objective function
        
        Returns
        -------
        numpy 2-d array where the columns corrospond to the objective functions
        and the rows corrospond to individuals for easy plotting
        '''
        data_list = self.getFront(gen_number, front_number, obj_functions, VFO)
        returned_data =np.array(np.zeros((len(data_list), len(obj_functions))))
        #This just data copying
        for i in range(0,len(data_list)):
            for j in range(1,len(data_list[i])):
                returned_data[i,j-1] = data_list[i][j]        

        return returned_data
        
        
        
    def windowFront(self, pts_list, window_list):
        '''
        Windows a front by removing any points that fall outside of a given
        bounding box.
        
        Parameters
        ----------
        pts_list:   list-of-list of points, with no leading ID number. Assumed
                    to come from something like reduceFront.
                    
        window_list:    list of 2-entry lists, where each 2-entry list 
                        corrosponds to an objective, and the 2-entry list
                        gives the minimum (bin[0]) and maximum (bin[1])
                        values for the list
        '''
        ret_list = []
        
        for pt in pts_list:
            #If point fall within window
            okpt = True
            #Zip together points with min-max so we can handle any number
            #of objectives
            for combined  in zip(pt, window_list):
                if (combined[0] < combined[1][0]):
                    okpt = False
                if (combined[0] > combined[1][1]):
                    okpt = False
            #If OK, then add
            if okpt:
                ret_list.append(pt)
        
        return ret_list
    
    def _UtilAppendFront(self, numpy_data, returned_front_data, scale_factor):
        #Moves data from a list of tuples into Numpy arrays for plotting
        #While stripping leading ID term
        #If numpy_data = none, a new numpy array will be created to be 
        #returned
        
        #Make up a numpy object 
        num_obj = len(returned_front_data[0])-1
        num_len = len(returned_front_data)
        
        new_data = np.empty((num_len, num_obj))
        
        for i in xrange(0, num_len):
            for j in xrange(0, num_obj):
                if returned_front_data[i][j + 1] != 'Infeasible':
                    new_data[i][j] = returned_front_data[i][j + 1]*scale_factor[j]
                    
        if (numpy_data == None):
            ret_data = new_data
        else:
            ret_data = np.vstack((numpy_data, new_data))
        
        return ret_data
        
            
    
    def GenPlot2D(self, gen_num, obj_functions, scale_factor, filename):
        '''Builds a single png of 2 or 3 objective function results
        from a single generation in the database  
        
        
        Parameters
        ----------
        gen_num:  Number of the generation to plot
    
        obj_functions: List-like listing of objective function IDs to pull
                       corrosponding objective function for individuals.
                       Will error if not either 2 or 3 long
        
        filename: filename for the final converted GIF
        
        Returns
        -------
        No return value, will create filename.png
        
        '''        
        #Check to see if num_obj is reasonable
        num_obj = len(obj_functions)
        if (num_obj != 2):
            print('Error - GenPlot can only plot 2d plots')
            raise Exception('opt - post_process - GenPlot argument error')
        if (num_obj != len(scale_factor)):
            print('Error - GenPlot requires a scale factor (can be  1.0)')
            print('For all objective functions')
            raise Exception('opt - post_process - GenPlot argument error')
            
        #Do the basic plotting set-up
        plt.xlabel('Objective ' + str(obj_functions[0]))
        plt.ylabel('Objective ' + str(obj_functions[1]))
        plt.title('Generation ' + str(gen_num))
 
        data_flag = True
        current_front = 0
        reserve_data = None 
        points_style = ['ro', 'b*', 'g^', 'yd', 'ks', 'mp', 'c+'] 
        legend = []
        while (data_flag):
            #Try to pull data
            data = self.getFront(gen_num, current_front, obj_functions)
            if (len(data) > 0) and (current_front < 7):
                #If data, see if we have a point style ready - only 
                #first 7 fronts get distinct styles
                self.addFrontTo2DPlot(data,
                           format_string = points_style[current_front])
                legend.append('Front ' + str(current_front))
            elif (len(data) > 0) and (current_front >= 7):
                #Otherwise, add to reserve data to be plotted with a 
                #single plot in a bit
                if (reserve_data == None):
                    reserve_data = data
                    legend.append('Remaining Fronts')
                else:
                    np.vstack((reserve_data, data))
            else:
                #No data, done with loop
                data_flag = False
            #Move the counter
            current_front += 1
        #Out of loop, if we had more than 7 fronts, plot the rest
        if (reserve_data != None):
            self.addFrontTo2DPlot(reserve_data, 'ko')
        plt.legend(legend)
                           
        plt.savefig(filename + '.png', dpi=144)
        plt.close()
    
    
    def KrigingStats(self, generations, ident, tol= 1.e-4):
        '''Prints out Kriging model statistics for all generations in the 
        generation list. 
        
        Statistics include the number of points, number 
        of points updated by the Kriging model, min, max, average, and 
        standard deviation of the Kriging update by front.  Overall
        min/max also listed for all fronts and generations.
        
        Parameters
        ----------
        generations: list-like list of generation numbers
        
        ident: objective function ID to pull kriging model corrections from
        
        tol: Tolerance - deviation from 1.0 required to count as non-corrected
             point.  If fabs(correction - 1.0 ) <= tol point is not considered
             corrected. 
    
        Returns
        -------
        Tuple of lists, where each list has one entry for each generation
        in the generations input.  Lists are: maximum correction, minimum 
        correction, average correction, standard deviation of correction, 
        and fraction of entries with correction
        other than 1.0 to tol.
        '''
        #Make up return lists
        max_corrections = []
        min_corrections = []
        avg_corrections = []
        std_corrections = []
        fraction_corrections = []
        
        for gen in generations:
            sel_tuple = (gen, ident)
            self.cur.execute('select * from objfun where Generation= ? AND ' + \
                'ObjID = ?', sel_tuple)
            #Total number of points
            num_pts = 0.
            #Scan the array of corrections and keep all those points tol away 
            #from 1.0
            correction_array = []
            for row in self.cur:
                num_pts += 1.
                if math.fabs(row[4] - 1.0) > tol:
                    correction_array.append(row[4])
            #Fill in the statistics for this generation
            if len(correction_array) > 0:
                correction_array = np.array(correction_array)
                max_corrections.append(correction_array.max())
                min_corrections.append(correction_array.min())
                avg_corrections.append(np.mean(correction_array))
                std_corrections.append(np.std(correction_array, ddof=1))
                fraction_corrections.append(float(len(correction_array))/
                                            float(num_pts))
            else:
                max_corrections.append(0.)
                min_corrections.append(0.)
                avg_corrections.append(0.)
                std_corrections.append(0.)
                fraction_corrections.append(0.)                
        return (max_corrections, min_corrections, avg_corrections,
                std_corrections, fraction_corrections)
    

    def ParetoStats(self, generations, span_norm=None):
        '''Prints out tracking statistics for the evolving Pareto front
        
        Statistics include number of points in Pareto front, average crowding
        distance, and optionaly a span metric defined as:
            product((max-min)/span_norm) for all objective functions

        Parameters
        ----------
        generations: list-like list of generation numbers to evaluate
        
        span_norm: optional, normalization factors for computing span of each
                   objective function.  If not passed, span vector will be
                   returned as none. Assumed to apply to all objective functions
                   in order.

    
        Returns
        -------
        Tuple of lists, where each list has one entry for each generation
        in the generations input.  Lists are: number of points in Pareto front,
        average crowding distance, and span metric if passed
        '''            
        
        #Make up return lists
        numpts = []
        avg_crowd = []
        span = None

        #do the points counter and average crowd metric
        for gen in generations:
            sel_tuple = (gen, 0)
            self.cur.execute('select * from generation where GenNum= ? AND ' + \
                'Front = ?', sel_tuple)            
            num_pts = 0.
            num_pts_crwd = 0. #Different to allow removable of infinite values
            crowd_dist_sum = 0.
            for row in self.cur:
                num_pts += 1
                if (row[3] != float('Inf')):
                    num_pts_crwd += 1. 
                    crowd_dist_sum += row[3]
            numpts.append(num_pts)
            if num_pts > 0:
                avg_crowd.append(crowd_dist_sum/num_pts_crwd)
            else:
                avg_crowd.append(0)
        
        #If requested, do the span normalization
        if (span_norm != None):
            span = []
            for gen in generations:
                spanfactor = 1.
                IDlist = []
                sel_tuple = (gen, 0)
                self.cur.execute(
                    'select * from generation where GenNum= ? AND ' + \
                    'Front = ?', sel_tuple) 
                for row in self.cur:
                    IDlist.append(row[1])
                obj_no = 1
                for norm in span_norm:
                    obj_vals =[]
                    for indiv in IDlist:
                        sel_tuple = (gen, indiv, obj_no)
                        self.cur.execute(
                        'select * from objfun where Generation= ?' 
                        + ' AND IndID= ? AND ObjID= ?', sel_tuple)
                        for row in self.cur:
                            obj_vals.append(row[3])
                    spanfactor *= math.fabs((max(obj_vals) - min(obj_vals))/
                                                     norm)
                    obj_no += 1
                span.append(spanfactor)
            
        return (numpts, avg_crowd, span)
        
        
    def getKrigingData(self, model_ID):
        '''
        Returns the Kriging x,y, sita data that corrosponds to model_ID
        
        Parameters
        ----------
        model_ID   Scalar model ID number for the Kriging model in the database
        
        Returns
        -------
        3-valued tuple, [0] is a numpy matrix of the x-data, [1] is the y
        observations in numpy matrix column format, and [2] is a numpy array
        of the sita values
        '''
        #Load the Kriging data
        sel_tuple = (model_ID,)
        self.cur.execute('select * from mmdata where modelID= ?', sel_tuple)
        IDList = []  #List of individuals in the model
        DataList = []  #List of corrosponding Kriging values
        for row in self.cur:
            IDList.append(row[1])
            DataList.append(row[2])
        
        #Figure out how long the chromosome is
        sel_tuple = (IDList[0],)
        self.cur.execute('select * from chromosomes where IndID=?', sel_tuple)
        count = 0
        #Sqlite uses sequential retrivial, so need to manually count rows!
        for row in self.cur:
            count += 1
        num_x_vals = count
        
        #Set up a numpy matrix to hold the data
        krig_x_vals = np.matrix(np.zeros((len(IDList), num_x_vals)))
        krig_y_vals = np.transpose(np.matrix(DataList)) #Should be column
        
        #Go through and get the X Vaues
        count = 0
        for indiv in IDList:
            sel_tuple = (indiv,)
            self.cur.execute('select * from chromosomes where IndID=?', 
                             sel_tuple)
            for row in self.cur:
                krig_x_vals[count, row[1]-1] = row[2]
            count +=1
        
        #Get the final sita values 
        sel_tuple = (model_ID,)
        self.cur.execute('select * from mmparameters where modelID=?', 
                         sel_tuple)
        sita = np.array(np.zeros(num_x_vals))
        for row in self.cur:
            sita[row[1] - 1] = row[2]
        
        return (krig_x_vals, krig_y_vals, sita)



    def getKriging(self, model_ID):
        '''
        Reconstructs a Kriging model from the database and returns a reference
        to the object
        
        Parameters
        ----------
        model_ID    Scalar model ID number for the Kriging model in the database
        '''
        #Get the data
        (krig_x_vals, krig_y_vals, sita) = self.getKrigingData(model_ID)
               
        #Make up the Kriging model
        recon_model = kriging.Kriging(krig_x_vals, krig_y_vals, sita)
        recon_model.setFixed(sita)
        
        return recon_model
    
    def KrigingPlot(self, model_ID, var_pairs, filename, size=20):
        '''
        Builds 3-D surface plots of Krigning models and plots to file
        
        A Kriging model ID must be provided from the database, along with 
        a series of pair of parameters to form the x/y plots, with the Kriging
        responses as the Z-axis

        Parameters
        ----------
        model_ID    Scalar model ID number for the Kriging model in the database

        var_pairs   List of 2-D tuples, each tuple specifies a gene ID to be
                    used as an x or a y range.  Range will automatically go 
                    from the max to the min of the values found in the
                    database during the plotting all other variables will be
                    held at their mean values
        
        filename    Base filename for output.  Actaul filename will be
                    basename_x(id)_y(id).png
                    
        size        number of points along each plotting axis.  Total number of
                    points in plot is size^2
        
        Returns
        -------
        None
        '''
        
        #Get the data
        (krig_x_vals, krig_y_vals, sita) = self.getKrigingData(model_ID)
               
        #Make up the Kriging model
        recon_model = kriging.Kriging(krig_x_vals, krig_y_vals, sita)
        recon_model.setFixed(sita)
        
        #Go through and do the plotting
        
        #Find max, min, and mean of columns
        mean_vals = np.mean(krig_x_vals, axis=0)
        max_vals = np.max(krig_x_vals, axis=0)
        min_vals = np.min(krig_x_vals, axis=0)
        
        #Do each pair of variables to plot
        for var in var_pairs:
        #Figure out the ranges of the variables
            print var
            step_x = (max_vals[0,var[0] - 1] - min_vals[0,var[0] -1])/ \
                                (float(size) - 1.)
            step_y = (max_vals[0,var[1] - 1] - min_vals[0,var[1] -1])/ \
                                (float(size) - 1.)
            
            #Make up plotting arrays
            plot_x = np.array(np.zeros((size)))
            plot_y = np.array(np.zeros((size)))
            for numstep in range(0, size):
                plot_x[numstep] = min_vals[0, var[0] - 1] + step_x * numstep
                plot_y[numstep] = min_vals[0, var[1] - 1] + step_y * numstep
            
            #Make a mesh of points x,y repeating over range
            print plot_x
            print plot_y
            x,y=np.meshgrid(plot_x, plot_y)
                    
            #Calculate the Kriging model at the points - code based on Jiandao's
            #initial work for prads 
            z=np.array(np.zeros([np.shape(x)[0],np.shape(x)[1]]))
    
            #X-value point
            plot_vector = np.matrix(mean_vals)
            
            for i in range(np.shape(x)[0]):
                for j in range(np.shape(x)[1]):
                    plot_vector[0, var[0] - 1] = x[i,j]
                    plot_vector[0, var[1] - 1] = y[i,j]
                    
                    #Confusing syntax - Kriging returns a matrix of tuples
                    #value, mse of prediction.  Get just the initial value
                    z[i,j]=recon_model.predictor(plot_vector)[0][0,0]
                    
            #Make the plot     
            fig=plt.figure()
            ax=Axes3D(fig)
            ax.plot_surface(x,y,z,rstride=1,cstride=1,cmap=cm.jet)
            ax.set_xlabel('Variable' + str(var[0]))
            ax.set_ylabel('Variable' + str(var[1]))
            ax.set_zlabel('Kriging Correction Factor')
            
            #Make up a filename
            filename_current = filename + 'x' + str(var[0]) + '_y' + \
                               str(var[1]) + '.png'
            plt.savefig(filename_current, dpi=144)
            plt.close('all')
        
        return
    
    def ObjMovie(self, start_gen, stop_gen, obj_functions, scale_factor, 
                 filename, refResult = None, ptLabel= 'NSGAII-VFO', forced_axis=None):
        '''Builds an animated GIF of 2 or 3 objective function results
        from a database run.  
        
        The GIF will have constant axis scaling. ImageMagik convert
        utility is required to be accessible on the command line to 
        build the GIF. 3-objective function plots will have two moves built,
        one of a 3-D projection of the points, and one a 3-graph series of 2-D
        cuts on O1/O2 O2/O3 O1/O3 axis
        
        Parameters
        ---------
        start_gen:  Scalar generation number to start with
        
        stop_gen: scalar generation to stop at
    
        obj_functions: List-like listing of objective function IDs to pull
                       corrosponding objective function for individuals.
                       Will error if not either 2 or 3 long
        
        filename: filename for the final converted GIF
        
        refResult:  2-D only, optional list of tuples for a reference curve to 
                    be plotted on each image. Tuple must contain the following:
                    0-index - list-like f1 values to plot
                    1-index - list-like f2 valuet to plot
                    2-index - string - plotting command (e.g. 'b-')
                    3-index - string - label for reference curve
        
        ptLabel:    Optional point label for the data points if refResult
                    is provided.  Defaults to "NSGAII-VFO"
        
        Returns
        -------
        No return value, will create filename.gif as a the animated gif
        
        '''
        
        #Internal parameters
        delay = 20 #Delay in hundreths of second between frames
        
        #Check to see if num_obj is reasonable
        num_obj = len(obj_functions)
        if (num_obj < 2) or (num_obj > 3):
            print('Error - ObjMovie can only plot 2d or 3d plots')
            raise Exception('opt - post_process - ObjMovie argument error')
        if (num_obj != len(scale_factor)):
            print('Error - ObjMovie requires a scale factor (can be  1.0)')
            print('For all objective functions')
            raise Exception('opt - post_process - ObjMovie argument error')
        
        #Set up axis bounds        
        axis_min = np.empty((num_obj))
        axis_max = np.empty((num_obj))
        
        #Data for all the plot frames
        frame_data = []
        
        #Go through and get the data 
        for gen in xrange(start_gen, stop_gen + 1):
            curr_front = 0  #Start with the Pareto front
            data = self.getFront(gen, curr_front, obj_functions)
            this_frame_data = self._UtilAppendFront(None, data, scale_factor)
            #Load data for next front
            curr_front = 1
            data = self.getFront(gen, curr_front, obj_functions)
            #While there is data, append to our existing data set
            while (len(data) > 0):
                this_frame_data = self._UtilAppendFront(this_frame_data, 
                                                         data, scale_factor)
                curr_front = curr_front + 1
                data = self.getFront(gen, curr_front, obj_functions)
            
            #Now figure out if we have to update axis min/max
            if (gen == start_gen):
                #Then no data - so set all data
                for i in xrange(0, num_obj):
                    axis_min[i] = np.amin(this_frame_data[:,i])
                    axis_max[i] = np.amax(this_frame_data[:,i])
            else:
                #We already have some min/max data
                for i in xrange(0, num_obj):
                    curr_min = np.amin(this_frame_data[:,i])
                    curr_max = np.amax(this_frame_data[:,i])                
                    if (curr_min < axis_min[i]):
                        axis_min[i] = curr_min
                    if (curr_max > axis_max[i]):
                        axis_max[i] = curr_max
            
            #Store in the correct place
            frame_data.append(this_frame_data)
        
        #Now do the plots
        index = 0
        for frame in frame_data:
            temp_filename = filename + '_gen_' + str('%05d' % index) + '.png'
            if (num_obj == 2):
                if refResult is not None:
                    for refSet in refResult:
                        plt.plot(refSet[0], refSet[1], refSet[2],
                             label=refSet[3])
                    plt.plot(frame[:,0], frame[:,1], 'ro', label=ptLabel)
                    plt.legend()
                else:
                    plt.plot( frame[:,0], frame[:,1], 'ro')
                if not forced_axis:
                    plt.axis([axis_min[0], axis_max[0], axis_min[1], axis_max[1]])
                else:
                    plt.axis(forced_axis)
                plt.xlabel('Objective ' + str(obj_functions[0]))
                plt.ylabel('Objective ' + str(obj_functions[1]))
                plt.title('Generation ' + str(start_gen + index))
                plt.savefig(temp_filename, dpi=144)
                plt.close('all')
            elif (num_obj == 3):
                fig = plt.figure()
                ax = Axes3D(fig)
                ax.scatter(frame[:,0], frame[:,1], frame[:,2], color = 'r')
                ax.set_xlabel('Generation ' + str(start_gen + index) +
                ' Objective' + str(obj_functions[0]))
                ax.set_ylabel('Objective ' + str(obj_functions[1]))
                ax.set_zlabel('Objective ' + str(obj_functions[2]))
                ax.set_xlim3d(axis_min[0], axis_max[0])
                ax.set_ylim3d(axis_min[1], axis_max[1])                
                ax.set_zlim3d(axis_min[2], axis_max[2])
                plt.savefig(temp_filename, dpi=144)
                plt.close('all')
            
                ##Do some 2-D shoots as well
                plt.subplots_adjust(hspace = 0.4)
                plt.subplot(311)
                temp_filename2 = filename + '_gen2D_' + str('%05d' % index) + \
                          '.png'
                plt.plot( frame[:,0], frame[:,1], 'ro')
                plt.axis((axis_min[0], axis_max[0], axis_min[1], axis_max[1]))
                plt.xlabel('Objective ' + str(obj_functions[0]))
                plt.ylabel('Objective ' + str(obj_functions[1]))
                plt.title('Generation ' + str(start_gen + index))  
                
                plt.subplot(312)
                plt.plot( frame[:,1], frame[:,2], 'ro')
                plt.axis((axis_min[1], axis_max[1], axis_min[2], axis_max[2]))
                plt.xlabel('Objective ' + str(obj_functions[1]))
                plt.ylabel('Objective ' + str(obj_functions[2]))
                
                plt.subplot(313)
                plt.plot( frame[:,0], frame[:,2], 'ro')
                plt.axis((axis_min[0], axis_max[0], axis_min[2], axis_max[2]))
                plt.xlabel('Objective ' + str(obj_functions[0]))
                plt.ylabel('Objective ' + str(obj_functions[2]))                
                
                
                
                plt.savefig(temp_filename2, dpi=144)  
                plt.close('all')                 
            else:
                print('Coding error - post process movie file')
                print('This line should not be reachable')
            index = index + 1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
        
        #Make the movie
        subprocess.check_call('convert -delay '+ str(delay) +\
                              ' ' + filename + '_gen_*.png ' +filename +'.gif',
                              shell=True)
        if (num_obj == 3):
            subprocess.check_call('convert -delay '+ str(delay) +\
                              ' ' + filename + '_gen2D_*.png ' +filename +
                              '2D.gif',
                              shell=True)            
        
