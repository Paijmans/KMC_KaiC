#!/usr/bin/env python

import sys

import numpy as np
import scipy.io
from scipy.special import erfc
import matplotlib.pyplot as plt
from math import *


'''
Usage:

> ./histnow1D.py <filename> <#bins> <ColumnToShow>

data file should end with '.dat'.

'''

def load_data(filename):
    
    try:      
        data_file = open(filename)                   
    except:
       print "No data file named %s" %sys.argv[i]
       sys.exit() 

    # Read one or more lines.
    columnnbr = len(data_file.readline().split())
    #data = np.array([float(x) for x in data_file.read().split()], np.float)   
    data = np.loadtxt(filename, delimiter='\t', skiprows=0, usecols=range(0,columnnbr))
      
    data_file.close()

    return np.transpose(data), columnnbr

def make_histogram( dataarray, columnnbr, binwidth ):
    
    returnarray=[]   
    for data in dataarray:
      
      print "Fraction with KaiA: %s" %(float(sum(data > 0))/len(data))
            
      xmin=min(data)
      xmax=max(data)
      binnbr = int((xmax-xmin)/binwidth)
      Ntot = len(data)
      sys.stderr.write("xmin: %s, xmax: %s, bin#: %s, for %s samples\n" %(xmin, xmax, binnbr, Ntot))
      
      hist1, rx1 = np.histogram(data,bins=binnbr, range=(xmin,xmax), normed=False)                    
      rx1b = [(rx1[i]+rx1[i+1])/2 for i in range(binnbr)]
    
      hist1norm = np.divide(hist1, Ntot*(rx1b[1]-rx1b[0]))

      sys.stderr.write("Mean = %s, Std: %s\n" %(np.mean( data ), np.std( data )))
      
      returnarray.append([rx1b, hist1norm])

    return returnarray

#Write histogram data to file.
def write_to_file( data, inputfile_name ):

    for i in range(len(data)):

      output_file_name = inputfile_name.replace('.dat','.hist'+str(i))
      try:      
         output_file = open(output_file_name, "wt")                   
      except:
         print "Error in creating %s" %output_file_name
         sys.exit() 

      output_file.write("#x\thist\n")
      for j in range(len(data[i][0])):
        output_file.write("%s\t%s\n" %(data[i][0][j],data[i][1][j]))

      output_file.close()


def plot_hist(hist):
 
    plt.figure(1)
    plt.title("Histogram")
    plt.plot(hist[0],hist[1],'-',hist[0],hist[1],'.')     

    plt.show()

    
if __name__ == '__main__':

    if( len(sys.argv[1:]) == 0 ):
        print "Please provide input file."
        sys.exit()
       
    try:
      binwidth = float(sys.argv[2])
    except:
      binwidth = 1.0
      
    try:
      showplot = int(sys.argv[3])
    except:
      showplot = -1

    sys.stderr.write("Loading Filename: %s, binwidth: %s\n" %(sys.argv[1], binwidth))
    
    data, columnnbr = load_data( sys.argv[1] )
    datahistograms  = make_histogram( data, columnnbr, binwidth )
    write_to_file( datahistograms, sys.argv[1] )
    
    if( showplot > -1 ):
      plot_hist(datahistograms[showplot])
    



