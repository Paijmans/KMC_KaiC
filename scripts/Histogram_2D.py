#!/usr/bin/env python

import sys

import numpy as np
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
    header = data_file.readline().split()    
    data = array([float(x) for x in data_file.read().split()], np.float)   
    data_file.close()

    columnnbr = len( header )
    col1 = data[column1-1::columnnbr]
    col2 = data[column2-1::columnnbr]   

    print "Reading columns %s and %s from total column nbr %s" %(column1, column2, columnnbr)
    
    return col1, col2
    

def make_2Dhist(data1, data2 ):

    print "column1#: %s, column2#: %s" %(len(data1), len(data2))
    
    #Make histogram of period lengths.
    xedges = [-6.5,-5.5,-4.5,-3.5,-2.5,-1.5,-0.5,0.5,1.5,2.5,3.5,4.5,5.5,6.5]
    yedges = [-0.5,0.5,1.5,2.5,3.5,4.5,5.5,6.5]
    hist, rx1, ry1 = np.histogram2d(data1, data2, bins=(xedges,yedges), normed=True)       

    return hist


def plot_hist( hist ):
    rx2 = range(0,7,1)
    ry2 = range(0,7,1)

    loghist = hist
    maxval=loghist.max()
    minval=loghist.min()

    print loghist

    c = plt.contourf(rx2,ry2,loghist,linspace(minval,maxval,100))
    b = plt.colorbar(c, orientation='vertical')
    plt.winter()
    lx = plt.xlabel("x")
    ly = plt.ylabel("y")
    ax = plt.axis([0,6,0,6])
    plt.show()
    
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
    data1, data2  = load_data( sys.argv[1], column1, column2 )
    hist = make_2Dhist( data1, data2 )

    print hist

    #Write output to file.
    of = open('fluxhist.dat', 'w')    
    of.write("npT\tnpS\tSProb(npT,npS)\n")
    for row1 in range(0,13,1):
      for row2 in range(0,7,1):
        of.write("%s %s %s\n" %(row1-6, row2, hist[row1][row2]))
   
    #plot_hist( data1, data2 )

    #xlabel(r'X', size=20)
    #ylabel(r'$p(X)$', size=20)
   
    #show()



