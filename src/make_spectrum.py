#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
Project Details: 
    make_spectrum.py 
    Creates SDSS_xid.csv and specObj.csv when a particular RA and DEC are given
    To be used alongside Bubbleimg
    
    Link to Bubbleimg: "https://github.com/aileisun/bubbleimg"

Created on Mon Aug 14 15:58:47 2017

__author__      = "nnarenraju"
__copyright__   = "Copyright 2017, NameOfProject"
__credits__     = "nnarenraju"
__license__     = "Apache License 2.0"
__version__     = "1.0.1"
__maintainer__  = "nnarenraju"
__email__       = "nnarenraju@gmail.com"
__status__      = "inUsage"

Github Repository: "https://github.com/nnarenraju/ASIAA_SSP"

"""

import os
import glob
import shutil
import numpy as np
import astropy.table as at
import astropy.coordinates as ac

from astropy.io import ascii
from PyAstronomy import pyasl
from astropy import units as u
from astropy.coordinates import Angle

def _initialise(location):
        """Returns astropy table of given path to FITS/CSV file.""" 
        return at.Table.read(location)

def __SAS__(ra, dec, constraint, table):
    """Search Around Sky function"""
    catalog = ac.SkyCoord(ra=table['ra']*u.degree, dec=table['dec']*u.degree)    
    scalarc = ac.SkyCoord(ra*u.degree, dec*u.degree)    
    d2d = scalarc.separation(catalog)  
    catalogmsk = Angle(d2d).arcsec < constraint
    idxc = np.where(catalogmsk)[0]
    return idxc

def make_spectrum(ra, dec, save_dir, spec_table, dr14_data, constraint=5, overwrite=False):
    """
    Make spectrum
    """
    #Initialising Spectrum Table
    table = _initialise(spec_table)
    
    def _cluster(cluster_path, constraint=5):
        """Make cluster"""
        constraint = float(constraint)
        cluster_table = _initialise(cluster_path)
        cluster = __SAS__(ra, dec, constraint, cluster_table)
        return cluster
    
    def _get_AngDist(pos):
        """Returns angular distance between given ra, dec and compared ra, dec"""
        AngDist=(Angle(pyasl.getAngDist(table['ra'][pos], table['dec'][pos], ra, dec), u.degree))
        return AngDist
        
    
    #Get Parameters
    def _get_pos():
        """Returns RUN2D and plate vales of req record"""
        cluster = _cluster(spec_table, constraint=constraint)
        if len(cluster)>1:
            print "More than one instance found, choosing closest one"
        elif len(cluster)==0:
            print "Sample without sciencePrimary"
            return -1
        AngDist = map(_get_AngDist, cluster)
        pos_dr14 = cluster[AngDist.index(min(AngDist))]
        return pos_dr14
    
    def _get_plate():
        pos_dr14 = _get_pos()
        plate=str(table[pos_dr14]['#plate'])
        #Correct Plate
        if int(plate)/1000==0:
            plate='0'+plate
        return plate
    
    def _get_fiberid():
        pos_dr14 = _get_pos()
        fiber=str(table[pos_dr14]['fiberid'])
        #Correct fiberid
        if int(fiber)/10==0:
            fiber='000'+fiber
        elif int(fiber)/100==0:
			fiber='00'+fiber
        elif int(fiber)/1000==0:
			fiber='0'+fiber
        return fiber
        
    def _get_mjd():
        pos_dr14 = _get_pos()
        mjd=str(table[pos_dr14]['mjd'])
        return mjd
    
    def _get_RUN2D():
        pos_dr14 = _get_pos()
        if pos_dr14==-1:
            return -1
        return str(table[pos_dr14]['run2d'])
    
    def _get_record():
        pos_dr14 = _get_pos()
        return table[pos_dr14]
        
    def _return_paths(database):
        """Returns a path based on given iterable"""
        RUN2D=_get_RUN2D()
        if RUN2D==-1: 
            return -1
        spectrum_path=dr14_data+'dr14/'+database+'/spectro'+'/redux/'+\
                        RUN2D.replace("'", "")+'/spectra'+'/lite'
        return spectrum_path
    
    #MAIN CODE
    if not os.path.isfile(save_dir+'spec.fits'):
            
        #Input parameters
        DB=['eboss', 'sdss']
        paths = map(_return_paths, DB)
        if -1 in paths:
            return False
        required_paths = filter(os.path.exists, paths)
        print "DR14 counterpart of Object present in: ", required_paths
        
        for path in required_paths:
            plate=_get_plate()
            mjd=_get_mjd()
            fiber=_get_fiberid()
            store_name = glob.glob(str(path)+'/'+plate+'/'+'spec-'+plate+'-'+mjd+'-'+fiber+'.fits')
            shutil.copy(src=store_name[0], dst=save_dir)
            name = glob.glob(save_dir+'/*')
            os.rename(name[0], save_dir+'/spec.fits')
            
            
        if not os.path.isfile(save_dir+'/sdss_xid.csv'):
            ascii.write(_get_record(), save_dir+'/sdss_xid.csv', delimiter=",")
        if os.path.isfile(save_dir+'/spec.fits'):
            print "\n\n"
            return True
        else:
            print "\n\n"
            return False
    else:
        print "[sdssObj]Skipped"
        return True
