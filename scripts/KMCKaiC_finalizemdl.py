#!/usr/bin/env python

import sys

import numpy
import scipy.io
from scipy.special import erfc
from matplotlib.pylab import *
from math import *
'''
Usage:

> ./histnow.py <project_name>

model name should not contain any extention.

'''

def find_between( s, first, last ):
    try:
        start = s.index( first ) + len( first )
        end = s.index( last, start )
        return s[start:end]
    except ValueError:
        return ""
        
        
def find_between_r( s, first, last ):
    try:
        start = s.rindex( first ) + len( first )
        end = s.rindex( last, start )
        return s[start:end]
    except ValueError:
        pass
        
        
    
if __name__ == '__main__':

    if( len(sys.argv[1:]) == 0 ):
        print "Please provide model name."
        sys.exit()

    project_name = sys.argv[1]
    orig_mdl_file = project_name + '.par'
    
    try:      
        mdl_file = open(orig_mdl_file)
    except:
        print "No model file named %s" %orig_mdl_file
        sys.exit()
        
    temp_file = open("out.temp", "wt")
        
    for line in mdl_file:
        
        arithmitic = find_between(line,'[',']')
        toreplace = '[' + arithmitic + ']'
        
        if(arithmitic != ''):
            try:
                temp_file.write(line.replace(toreplace, str(eval(arithmitic)) ))
            except ValueError:
                print "No valid arithmitic: %s" %arithmitic
                sys.exit()
        else:
            temp_file.write( line )
        
    temp_file.close()      
    

'''
Fancy functions/ Greens functions:

1/sqrt(2*pi*sigma_sq)*exp( -x*x/(2*sigma_sq) )

'''

    #xlim(0.9, 2.2e2)
    #ylim(2e-6, 2e1)
    #xticks([1, 10, 100], ['1', '10', '100'], size=22)
    #yticks(size=18)
    #solline.set_label(r'theory')
    #legend(handlelen=0.02, pad=0.02,handletextsep=0.01, labelsep=0.001)
    #grid()
