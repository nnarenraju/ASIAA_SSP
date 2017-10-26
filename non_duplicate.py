#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
Project Details: 
    non_duplicate.py 
    Retains a single instance of a galaxy, if more than one is present 

Created on Sun Oct 15 22:54:58 2017

__author__      = "nnarenraju"
__copyright__   = "Copyright 2017, AGN_Classification"
__credits__     = "nnarenraju"
__license__     = "GPL"
__version__     = "1.0.1"
__maintainer__  = "nnarenraju"
__email__       = "nnarenraju@gmail.com"
__status__      = "inUsage"

Github Repository: "Add link here"

"""

import os
import astropy.table as at

from astropy.io import ascii
from astropy import units as u
from collections import defaultdict
from astropy.coordinates import SkyCoord


def _get_instance(cluster, table, **kwargs):
    """ Eliminates unrequired instances """
    priority = []
    for instance in cluster:
        priority.append(kwargs[table["Database"][instance]])
    return cluster[priority.index(min(priority))]
    
def _initialise(location):
    """ Returns an astropy table if possible """
    if not os.path.exists(location):
        raise NameError("File not present in given location")
    try: 
        return at.Table.read(location)
    except: 
        raise NameError("File cannot be opened as an astropy table")

def __SAS__(constraint, cluster_table, radec_table=False):
    """ Search Around Sky Function """
    if ~radec_table:
        radec_table=cluster_table
    radec_table = SkyCoord(ra=radec_table['RA']*u.degree, dec=radec_table['DEC']*u.degree)
    cluster_table = SkyCoord(ra=cluster_table['RA']*u.degree, dec=cluster_table['DEC']*u.degree)
    idxc1, idxc2, _, _ = cluster_table.search_around_sky(radec_table, constraint)
    return idxc1, idxc2    
    
def _cluster(cluster_path, radec_path=False, constraint=5):
        """ 
        Clusters the spectra obtained with the 'n' arcsec constraint.
        """
        def _list_duplicates(idxc1, idxc2):
            """ Returns positions of similar items """
            clusters = []
            tally = defaultdict(list)
            for i,item in enumerate(idxc1):
                tally[item].append(i)
            for key,locs in tally.items():
                clusters.append([idxc2[loc] for loc in locs])
            clusters = [sorted(cluster) for cluster in clusters]
            return [list(x) for x in set(tuple(x) for x in clusters)]
        
        constraint=float(constraint)*u.arcsec
        cluster_table = _initialise(cluster_path)
        if ~radec_path:
            radec_table = False
        else:
            radec_table = _initialise(radec_path)
        
        # Get similarity using Search Around Sky
        idxc1, idxc2 = __SAS__(constraint, cluster_table=cluster_table, radec_table=radec_table)
        clusters = _list_duplicates(idxc1, idxc2)
        
        if len(clusters):
            return sorted(clusters)

if __name__ == "__main__":
    # Input table to be non-duplicated
    table_path = 'Unified_database.csv'
    table = _initialise(table_path)
    # Cluster the instances
    clusters = _cluster(cluster_path=table_path)
    
    # Non-duplicated instances based on given priority
    non_dupe = []
    for cluster in clusters:
        non_dupe.append(_get_instance(cluster, table, Mullaney=1, Zakamska=2, Reyes=3, Yuan=4))

    # Get unrequired rows and remove them
    _list = [i for i in range(len(table)) if i not in non_dupe]
    table.remove_rows(_list)

    # Write required data into a CSV file
    ascii.write(table, "non_duplicated.csv", delimiter=",", overwrite=True)
    