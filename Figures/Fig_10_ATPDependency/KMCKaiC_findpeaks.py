#!/usr/bin/env python

#To run type:
# > KMCKaiC_findperiods <full_filename> <'1' for plot>

import sys

import numpy as np
import scipy.io
from scipy.special import erfc
from matplotlib.pylab import *
from math import *
from Pypeakfinder import *

if __name__ == '__main__':

    if( len(sys.argv[1:]) == 0 ):
        print "Please provide KMCKaiC timetrace file"
        print "USE: >KMCKaiC_findperiods <full_filename> <'1' for plot>"
        sys.exit()

    datafile_name = sys.argv[1]
    showplot = -1
                    
    try:
        showplot = int(sys.argv[2])
    except:
        pass
    
        
    outputmin_name = datafile_name.replace('.0.','.') + '.minima'
    outputmax_name = datafile_name.replace('.0.','.') + '.maxima'
    
    try:      
        data_file = open(datafile_name,'r')                 
    except:
        print "No data file named %s" %datafile_name
        sys.exit()
        
    try:      
        outputmax_file = open(outputmax_name,"wt")              
        outputmin_file = open(outputmin_name,"wt")
    except:
        print "Something went wrong in opening min/max output file named %s / %s" %outputmax_name, outputmin_name
        sys.exit()
    
    header = data_file.readline().split()
    data = array([float(x) for x in data_file.read().split()], np.float)
    #data = np.loadtxt(datafile_name, skiprows=1, delimiter="\t")
      
    data_file.close()
    
    column_nbr = len(header)
   
    times = data[0::column_nbr]
    phosp = data[1::column_nbr]
    activ = data[22::column_nbr]

    dt = times[2]-times[1]

    #Lookahead window is now 6h (6/0.01=600).
    lookahead_val = 6/dt
      
    [maxima, minima] = peakdetect(y_axis = phosp, x_axis = times, lookahead = lookahead_val, delta = 0)
    maxima = np.array(maxima)
    minima = np.array(minima)    
    
    
    #Writing maxima/ minima to files
    for i in range(min(len(maxima),len(minima))):
        outputmax_file.write(str(str(maxima[i][0]) + '\t' + str(maxima[i][1]) + '\n'))
        outputmin_file.write(str(str(minima[i][0]) + '\t' + str(minima[i][1]) + '\n'))
    
    outputmax_file.close()
    outputmin_file.close()    
    
    #Data analysis
    PtPtimes=[]
    Amplitudes=[]
    for i in range(len(maxima)-1):
        PtPtimes.append(maxima[i+1][0]-maxima[i][0])
        Amplitudes.append(0.5*(maxima[i][1]-minima[i][1]))

    tphostimes=[]
    tdephtimes=[]        
    if( maxima[0][0] < minima[0][0] ):
        for i in range( min( len(maxima), len(minima) ) - 1 ):
            tphostimes.append( maxima[i+1][0] - minima[i][0] )
            tdephtimes.append( minima[i][0] - maxima[i][0] )
    else:
        for i in range( min( len(maxima), len(minima) ) - 1 ):
            tphostimes.append( maxima[i][0] - minima[i][0] )
            tdephtimes.append( minima[i+1][0] - maxima[i][0] )
        
    MeanPeriod = np.mean(PtPtimes)
    StdPeriod = sqrt(np.mean(np.square(PtPtimes)) - MeanPeriod**2)
    
    MeanAmplitudes = np.mean(Amplitudes)
    StdAmplitudes = sqrt(np.mean(np.square(Amplitudes)) - MeanAmplitudes**2)
    
    Meantphos=np.mean(tphostimes)
    Meantdeph=np.mean(tdephtimes)
       
    #for i in range( min( len(tphostimes), len(tdephtimes) ) - 1 ):
    #    print tphostimes[i], tdephtimes[i]
        
    MeanConc = np.mean(phosp)
    MeanMin = np.mean(minima[:,1])
    MeanMax = np.mean(maxima[:,1])
    
    print MeanPeriod, StdPeriod, MeanAmplitudes, StdAmplitudes, MeanMin, MeanMax, MeanConc, Meantphos, Meantdeph
    
    if showplot==1:
    
        #For plotting:
        maxtimes=[]
        maxconcs=[]
        plot(times,phosp,'-')
        plot(maxima[:,0],maxima[:,1],'o')
        plot(minima[:,0],minima[:,1],'o')        
        plt.show()
        
        print str(len(PtPtimes))
        hist, bin_edges = np.histogram(PtPtimes, density=True,bins=len(PtPtimes)/10)
        plot(hist)
        plt.show()


