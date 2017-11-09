#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
Project Details: 
    HSC_filter.py 
    Plotting the HSC - Subaru Filter Transmission Function

Created on Mon Oct 16 02:36:47 2017

__author__      = "nnarenraju"
__copyright__   = "Copyright 2017, AGN_Classification"
__credits__     = "nnarenraju"
__license__     = "Apache License 2.0"
__version__     = "1.0.1"
__maintainer__  = "nnarenraju"
__email__       = "nnarenraju@gmail.com"
__status__      = "inUsage"

Github Repository: "https://github.com/nnarenraju/ASIAA_SSP"

"""

import os
import numpy as np
import astropy.table as at
import matplotlib.pyplot as plt


def _initialise(location):
    """ Returns an astropy table if possible """
    if not os.path.exists(location):
        raise NameError("File not present in given location")
    try: 
        return at.Table.read(location)
    except: 
        raise NameError("File cannot be opened as an astropy table")
        
def _plot_func(table, band, marker, heuristic_param=0.30, pos=1):
    """ Plot transfer function """
    param = heuristic_param
    plt.plot(table['ws'], table['trans']/max(table['trans']), label = "{}".format(band))
    for x, y in zip(table['ws'], table['trans']/max(table['trans'])):
        if pos==2 and x < np.median(table['ws']):
            continue
        elif abs(y-param*max(table['trans']/max(table['trans'])))<0.1:
            plt.plot(x, y, marker=marker, label='{0}% of {1}'.format(param, band))
            print "{0} % of {1}".format(param, band),(x,y)
            break
        
if __name__ == "__main__":
    # initialise all required data
    parent_dir = ""
    filenames = ["g.csv", "r.csv", "i.csv", "z.csv", "y.csv"]
    filenames = [parent_dir+filename for filename in filenames]
    _g, _r, _i, _z, _y = map(_initialise, filenames)
    
    # Plot individial band transfer function
    plt.figure(figsize = (30,10))
    # In the order = g, r, i ,z, y and y_end
    _plot_func(_g, "g", "x", heuristic_param=0.30)
    _plot_func(_r, "r", "o", heuristic_param=0.76)
    _plot_func(_i, "i", "*", heuristic_param=0.56)
    _plot_func(_z, "z", "s", heuristic_param=0.32)
    _plot_func(_y, "_y", "D", heuristic_param=0.50)
    _plot_func(_y, "y_", "^", heuristic_param=0.30, pos=2)    
    
    plt.xlabel("Frequency")
    plt.ylabel("Amplitude")
    plt.grid(True)
    plt.legend()
    plt.title("HSC Filter Transmission Function")
    plt.show()
