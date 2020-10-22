#!/usr/bin/env python

import sys

import numpy as np
import math
import scipy.io
from scipy.special import erfc
from matplotlib.pylab import *
from math import *
'''
Usage:

> ./histnow.py <file1> <file2> ... <#bins>

data file should end with '.dat'.

'''

def load_data(filename, column1, column2):
    
    try:      
        data_file = open(filename)                   
    except:
       print "No data file named %s" %sys.argv[i]
       sys.exit() 

    # Read one or more lines.
    #header = data_file.readline().split()    
    data = array([float(x) for x in data_file.read().split()], np.float)   
    data_file.close()

    columnnbr = 7 #len( header )
    col1 = data[0::columnnbr]
    col2 = data[1::columnnbr]
    col3 = data[2::columnnbr]   
    col4 = data[3::columnnbr]              
    col5 = data[4::columnnbr]   
    col6 = data[5::columnnbr]              
    col7 = data[6::columnnbr]              

    print "Reading columns %s and %s from total column nbr %s" %(column1, column2, columnnbr)
    
    return [col1, col2, col3, col4, col5, col6, col7]
    

def make_2Dhist( columns ):
    s1 = columns[0]
    s2 = columns[1]
 
    part_res=1.
    state1_nbr = int(max(s1)*part_res+1)
    state2_nbr = int(max(s2)*part_res+1)
    state1_max = max(s1)
    state2_max = max(s2)

    state_histogram = np.reshape( columns[2], (state1_nbr, state2_nbr) )
    fluxmap_s1p1 = np.reshape( columns[3], (state1_nbr, state2_nbr) )
    fluxmap_s1m1 = np.reshape( columns[4], (state1_nbr, state2_nbr) )
    fluxmap_s2p1 = np.reshape( columns[5], (state1_nbr, state2_nbr) )
    fluxmap_s2m1 = np.reshape( columns[6], (state1_nbr, state2_nbr) )  
    
    print "\n" 
    # Create output file
    off = open('fluxhist.dat', 'w')    
    off.write("npT\tnpS\tProb(npT,npS)\tfluxT\tfluxS\n")        
  
    dfluxmap_s1=np.zeros((state1_nbr,state2_nbr), dtype=np.float64)
    dfluxmap_s2=np.zeros((state1_nbr,state2_nbr), dtype=np.float64)    
   
    for i1 in range(0,state1_nbr,1):
      for i2 in range(0,state2_nbr,1):
        if( i1 < state1_max and i2 < state2_max ):
          dfluxmap_s1[i1][i2] = fluxmap_s1p1[i1][i2] - fluxmap_s1m1[i1+1][i2]
          dfluxmap_s2[i1][i2] = fluxmap_s2p1[i1][i2] - fluxmap_s2m1[i1][i2+1]
            
        elif( i1 == state1_max and i2 < state2_max ):         
          dfluxmap_s2[i1][i2] = fluxmap_s2p1[i1][i2] - fluxmap_s2m1[i1][i2+1]  
                
        elif( i1 < state1_max and i2 == state2_max ):          
          dfluxmap_s1[i1][i2] = fluxmap_s1p1[i1][i2] - fluxmap_s1m1[i1+1][i2]
                  
        totflux_s1 = 0.5*(dfluxmap_s1[i1][i2] + dfluxmap_s1[i1-1][i2])
        totflux_s2 = 0.5*(dfluxmap_s2[i1][i2] + dfluxmap_s2[i1][i2-1])                              

        #print i1/part_res, i2/part_res, hist_val, dfluxmap_s1[i1][i2], dfluxmap_s1[i1-1][i2]

        print i1/part_res, i2/part_res, state_histogram[i1][i2], totflux_s1, totflux_s2             
        off.write("%s\t%s\t%s\t%s\t%s\n" %(i1/part_res, i2/part_res, state_histogram[i1][i2], totflux_s1, totflux_s2))
        
            
if __name__ == '__main__':

    if( len(sys.argv[1:]) == 0 ):
        print "Please provide input file."
        sys.exit()

    try: 
        column1 = int(sys.argv[len(sys.argv[1:])-1])        
        column2 = int(sys.argv[len(sys.argv[1:])])
    except:
        column1 = 1; column2 = 2;
        
    print "Loading Filename: %s" %(sys.argv[1])
    columns = load_data( sys.argv[1], column1, column2 )
    make_2Dhist( columns )
   

