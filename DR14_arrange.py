#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
Project Details: 
    DR14_arrange.py 
    Arranges downloaded data from SDSS-DR14 into structure similar to 'bubbleimg'
    
    Link to Bubbleimg: "https://github.com/aileisun/bubbleimg"

Created on Fri Aug 11 11:38:47 2017

__author__      = "nnarenraju"
__copyright__   = "Copyright 2017, AGN_Classification"
__credits__     = "nnarenraju"
__license__     = "Apache License 2.0"
__version__     = "1.0.1"
__maintainer__  = "nnarenraju"
__email__       = "nnarenraju@gmail.com"
__status__      = "inComplete"

Github Repository: "https://github.com/nnarenraju/ASIAA_SSP"

"""

import os
import glob
import functools
import astropy.table as at
import astropy.coordinates as ac

from functools import reduce
from astropy.io import ascii
from PyAstronomy import pyasl
from astropy import units as u
from astropy.coordinates import Angle
from distutils.dir_util import copy_tree

class Make_readable(object):
    
    """ 
    Re-organises the DR14 default directory structure into a structure similar
    to that of bubbleimg.
    Groups together the spectrum files that are within 5 arcsec of the instance
    If more than one spectrum is available for a given instance, the closest 
    spectrum wrt angular distance is chosen and moved into the good directory. 
    If there exists no spectrum within 5 arcsec of a given instance, the 
    instance will be put into the except directory. 
    """
    
    def __init__(self):
        """Current status of the Object."""
        
        self._status = False
        self._make_combined_table = False
        self._clustering = False
        self._make_spectrum = False
        self._make_table = False
    
    def _initialise(self, location):
        """Returns astropy table of given path to FITS/CSV/txt file.""" 
        if "txt" in location: 
            return ascii.read(location)
        else: 
            return at.Table.read(location)

    def _combine(self, table, append_table):
        """Returns the concatenated astropy tables"""
        if type(table)!=at.Table or type(append_table)!=at.Table:
                raise TypeError('Input table is not an astropy table')
        return at.join(table, append_table, join_type='outer')
    
    def combine_tables(self, filename, dir_path=os.getcwd(), name='combined'):
        """Combine all tables that contain the filename as part of its name."""
        location = glob.glob(dir_path+"*"+filename+"*")
        if location == None:
            raise NameError('Table not located in the given path')
        else:
            tables = map(self._initialise, location)
            #Sanity check 0
            if type(tables[0])!= at.Table:
                raise TypeError('Astropy table could not be created')
            combined_table = reduce(self._combine, tables)
            if len(combined_table):
                self._make_combined_table = True
                try: 
                    combined_table.writeto(name+'.fits')
                except:
                    txt_file=open('download_rsync.txt', 'w')
                    c = combined_table
                    _list = [i for i in c[c.colnames[0]]]
                    txt_file.write("%s\n" %_list)
                return glob.glob(dir_path + name +'.fits')
                
    def __SAS__(constraint, cluster_table, radec_table=False):
        """Search Around Sky function"""
        def __SkyCoord__(table):
            table = ac.SkyCoord(ra = table['ra']*u.degree,
                                        dec = table['dec']*u.degree)    
            return table
        if ~radec_table:
            radec_table=cluster_table
        radec_table = __SkyCoord__(radec_table)
        cluster_table = __SkyCoord__(cluster_table)
        idxc1,idxc2,_,_=cluster_table.search_around_sky(radec_table,constraint)
        return idxc1, idxc2
        
    def _cluster(self, cluster_path, radec_path=False, constraint=5):
        """ 
        Clusters the spectra obtained with the 'n' arcsec constraint.
        table_path = Path of congregated spectra table list
            
        NOTE: Use make_spectrum_list to congregate Spectra Tables.
        """
        constraint=float(constraint)*u.arcsec
        cluster_table = self._initialise(cluster_path)
        if ~radec_path:
            radec_table = False
        else:
            radec_table = self._initialise(radec_path)
        idxc1, idxc2 = self.__SAS__(constraint, cluster_table=cluster_table, 
                                                     radec_table=radec_table)
        
        cluster = list(set(idxc1))
        for x in range(len(cluster)):
            index=[i for i,j in enumerate(idxc1) if j == cluster[x]]
            cluster[x].append([idxc2[i] for i in index])
            
        if len(cluster):
            self._clustering=True
            return cluster    

    def make_spectrum(self, table_path, input_table_path, dir_path=os.getcwd()):
        """
        Using the clusters formed using the _cluster definition this definition
        compares the input table with the clustered spectra table and moves the
        instances without a spectrum to except.
        """
        
        def _return_paths(database, RUN2D):
            """Returns a path based on given iterable"""
            spectrum_path=dir_path+'/dr14/'+database+'/spectro'+'/redux/'+ \
                          RUN2D+'/spectra'+'/lite/'
            return spectrum_path
        
        def __mkdir__(dir_path, dir_name):
            """Creates requested directory and changes current directory"""
            if not os.path.exists(dir_path+dir_name):
                os.makedirs(dir_path+dir_name)
                os.system("cd "+dir_path+dir_name)
        
        def __sanity__(path):
            """Check consistency for number of files/directories"""
            #Sanity Check 0
            files = folders = 0
            for _, dirnames, filenames in os.walk(path):
                files += len(filenames)
                folders += len(dirnames)
            return files, folders
        
        def __mkobj__(element, good=False):
            """Make directory for object eg., SDSSJXXXX:XXXX/CSV"""
            #Object_name should be a string object
            if not os.path.exists(object_name): 
                os.makedirs(object_name)
                if good:
                    return __mkgood__(element, object_name)
                else:
                    try: 
                        ascii.write(element, 'sdss_xid.csv')
                        return True
                    except: 
                        return False
                    
        def __mkgood__(element, object_name):
            """Move folder from DR14 using plate ID"""
            pos_cluster = clusters[clusters.index(element)][1]
            pos_radec = clusters[clusters.index(element)][0]
            AngDist = []
            for pos in pos_cluster:
                AngDist.append(Angle(pyasl.getAngDist(table['ra'][pos],
                                                      table['dec'][pos],
                                                      table['ra'][pos_radec],
                                                      table['dec'][pos_radec]),
                                                      u.degree))
            pos_dr14 = pos_cluster[AngDist.index(min(AngDist))]
            plate = table[pos_dr14]['plate']
            for path in required_paths:
                try:
                    functools.partial(copy_tree,dst=object_name),path+'/'+plate
                    stored_path = object_name+'/'+plate
                except:
                    continue
            if os.path.exists(stored_path):
                return True
            else:
                return False
            
        
        def _move_to_except(dir_path):
            """Handles the except directory sub-section"""
            __mkdir__(dir_path, dir_name='/except')
            list_except = input_table[:]
            not_req = map(lambda t: t[0], clusters)
            ascii.write(list_except.remove_rows(not_req), 'list_except.csv')
            status_mkobj = map(__mkobj__, list_except)
            if all(status_mkobj):
                print 'Directories have been succesfully created'
            else:
                raise ValueError('sdss_xid.csv not written for all objects')
            files, folders = __sanity__()
            if files==1 and folders==len(list_except):
                print "list_except.csv and req directories succesfully written"
            else:
                raise ValueError("Inconsistent number of files/directory")
            
        def _move_to_good(dir_path):
            """Handles the good directory sub-section"""
            __mkdir__(dir_path, dir_name='/good')
            req = map(lambda t: t[0], clusters)
            list_good = input_table[:]
            list_good = filter(lambda t: t.index in req, list_good)
            ascii.write(list_good, 'list_good.csv')
            status_mkobj = map(functools.partial(__mkobj__, good=True), req)
            if all(status_mkobj):
                print 'Directories have been succesfully created'
            else:
                raise ValueError('sdss_xid.csv not written for all objects')
            files, folders = __sanity__()
            if files==1 and folders==len(list_good):
                print "list_good.csv and req. directories succesfully written"
            else:
                raise ValueError("Inconsistent number of files/directory")
            
        #Downloading Spectra from SDSS DR14 Webpage using rsync
        command='rsync -avzL --files-from=download_rsync.txt \
                 rsync://data.sdss.org/dr14  dr14'
        os.system('cd '+dir_path+'/download_rsync.txt')
        os.system(command)
        
        #Input parameters
        DB=['eboss', 'sdss']
        RUN2D=['v5_10_0', '26', '103', '104']
        path=map(map(functools.partial(_return_paths, database=DB), RUN2D), DB)
        paths = map(lambda l:[item for sublist in l for item in sublist], path)
        required_paths = filter(os.path.exists, paths)
        clusters = self._cluster(cluster_path=table_path, 
                                 radec_path=input_table_path, 
                                 constraint=2)
        table = self._initialise(table_path)
        input_table = self._initialise(input_table_path)
        
        #Data Manipulation
        ascii.write(input_table, dir_path+'list.csv')
        _move_to_except(dir_path=dir_path)
        _move_to_good(dir_path=dir_path)
        
    def __str__(self):
        """Returns current status of object"""
        s1='Overall Status: {0}\n'.format(str(self._status))
        s2='Combined Spec Tables: {0}\n'.format(str(self._make_combined_table))
        s3='Clustering Process: {0}\n'.format(str(self._clustering))
        s4='Make Spectrum: {0}\n'.format(str(self._make_spectrum))
        s5='Re-organise Files/directories: {0}\n'.format(str(self._make_table))
        return s1, s2, s3, s4, s5
