#!/usr/bin/env python

#To run type:
# > Gillespie_findperiods <full_filename> <Td> <beta_c> <'1' for plot>

import sys

import numpy
import scipy.io
from scipy.special import erfc
from matplotlib.pylab import *
from math import *
from Pypeakfinder import *

if __name__ == '__main__':

    if( len(sys.argv[1:]) == 0 ):
        print "Please provide Gillespie timetrace file"
        sys.exit()

    datafile_name = sys.argv[1]
    Td     = -1
    beta_c = -1
    showplot = -1

    try:
        Td = float(sys.argv[2])
    except:
        pass 
        
    try:
        beta_c = float(sys.argv[3])
    except:
        pass
                      
    try:
        showplot = int(sys.argv[4])
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
    data = array([float(x) for x in data_file.read().split()], numpy.float)
    #d = numpy.loadtxt(datafile_name, skiprows=1, delimiter="\t")
      
    data_file.close()
   
    times = data[0::3]
    concs = data[1::3]
    volum = data[2::3]

    #Lookahead window is now 6h (6/0.01=600).
    [maxima, minima] = peakdetect(y_axis = concs, x_axis = times, lookahead = 600, delta = 0)
    
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
        
    AvgPeriod = numpy.mean(PtPtimes)
    StdPeriod = sqrt(numpy.mean(numpy.square(PtPtimes)) - AvgPeriod**2)
    
    AvgAmplitudes = numpy.mean(Amplitudes)
    StdAmplitudes = sqrt(numpy.mean(numpy.square(Amplitudes)) - AvgAmplitudes**2)
    
    MeanConc = numpy.mean(concs)
    MeanPrNbr = numpy.mean(numpy.multiply(concs, volum))
    
    print Td, beta_c, AvgPeriod, StdPeriod, AvgAmplitudes, StdAmplitudes, MeanPrNbr, MeanConc
    
    if showplot==1:
    
        #For plotting:
        maxtimes=[]
        maxconcs=[]
        for maximum in maxima:
            maxtimes.append(maximum[0])
            maxconcs.append(maximum[1])
        plot(times,concs,'-')
        plot(maxtimes,maxconcs,'o')
        plt.show()

	print len(PtPtimes)
	hist, bin_edges = np.histogram(PtPtimes, density=True,bins=50)
	plot(hist)
	plt.show()


