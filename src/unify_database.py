#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
Project Details:  
    unify_database.py 
    Combines the 4 required databases with required columns and appropriate naming 
    Creates a CSV and FITS file of the unified data

Created on Sun Oct 15 17:44:43 2017

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

import csv
import numpy as np
import astropy.table as at
import matplotlib.pyplot as plt
import warnings

from astropy.table import vstack
from astropy.io import ascii

"""
DATABASE INITIALISATION
"""
# Reading FITS file
parent_dir = ""
filenames = ["mullaney.fits", "zakamska.fit", "reyes.csv", "yuan.fits"]
file_locations = [parent_dir+filename for filename in filenames]
mullaney,zakamska,reyes,yuan = map(at.Table.read, file_locations)

"""
CALCULATING log10(OIII)
"""

# Yuan Database: OIII Luminosity present in req. form
""" Nothing to do """

# Reyes database: Convert LOIII from Sol Lum to erg/s
reyes_logOIII = reyes['LO3g'] + 33.58275

# Zakamska database: Convert LOIII from Sol Lum to erg/s
# Present in logarithmic form - need not convert
zakamska_logOIII = zakamska['L_OIII_'] + 33.58275

# Mullaney database: LOIII = OIII_5007_LUM + OIII_5007B_LUM
# Select AGN_TYPE = 2 (samples of concern)
mullaney = mullaney[mullaney['AGN_TYPE']==2]
LOIII = mullaney['OIII_5007_LUM'] + mullaney['OIII_5007B_LUM']
mullaney_logOIII = np.log10(LOIII)



"""
PLOTTING TYPE-II AGN DATABASE
"""

# Concatenating log10(OIII) column to corresponding databases
reyes['logOIII'] = reyes_logOIII
mullaney['logOIII'] = mullaney_logOIII
zakamska['logOIII'] = zakamska_logOIII

""" Plotting Redshift VS LOIII """

plt.figure(figsize=(12,8))
# Zakamska database
plot_Z = plt.scatter(zakamska['z'], zakamska['logOIII'], s = 1, label="Zakamska") 

# Mullaney database
plot_M = plt.scatter(mullaney['Z'], mullaney['logOIII'], s=1, color = 'green', label="Mullaney")

# Reyes database
plot_R = plt.scatter(reyes['z'], reyes['logOIII'], s = 1, color = 'red', label="Reyes")

# Yuan database
plot_Y = plt.scatter(yuan['z'], yuan['LOIII'], s = 1, color = 'orange', label="Yuan")

plt.xlabel("Stuff")
plt.ylabel("Amplitude")
plt.grid(True)
plt.legend()
plt.title("Unified Database")
plt.show()


"""
DATABASE UNIFICATION - Add Required fields
"""

def _cleanup(table, **kwargs):
    """ Concatenate required columns with appropriate naming """
    def _rename(new_name):
        """ Renames required column """
        # key=new_name, kwargs[key]=old_name
        if new_name!="Database":
            table.rename_column(kwargs[new_name], new_name)
    
    keep_cols = []
    # Rename required columns
    map(_rename, kwargs)
    map(keep_cols.append, kwargs)

    # Add Database column to table
    if "Database" in kwargs:
        table["Database"]=[kwargs["Database"]]*len(table)
    else:
        warnings.warn("Database column will not be added, readability might reduce", UserWarning)
    
    # Keep required columns
    table.keep_columns(keep_cols)
    return table

# Clean-up all databases
# Give appropriate old column names
yuan=_cleanup(yuan, Database="Yuan", RA="ra", DEC="dec", z="z", logOIII="LOIII")
mullaney=_cleanup(mullaney, Database="Mullaney", RA="RA", DEC="DEC", z="Z", logOIII="logOIII")
zakamska=_cleanup(zakamska, Database="Zakamska", RA="_RA", DEC="_DE", z="z", logOIII="logOIII")
reyes=_cleanup(reyes, Database="Reyes", RA="RAJ2000", DEC="DEJ2000", z="z", logOIII="logOIII")
    
unified_database = vstack([mullaney, zakamska, reyes, yuan])


# Writing unified database into a CSV and FITS file
with open('Unified_database.csv', 'wb') as f:
    writer = csv.writer(f)
    writer.writerows(unified_database)
    
ascii.write(unified_database, "Unified_database.fits", overwrite=True)
