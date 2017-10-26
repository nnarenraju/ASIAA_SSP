#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
Project Details: 
    HSC_convolve.py 
    Performning convolution operation to obtain AB Magnitude of Spectrum

Created on Mon Oct 16 03:27:52 2017

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
import numpy as np

from scipy import interpolate


def _initialise(location):
    """ Returns an astropy table if possible """
    if not os.path.exists(location):
        raise NameError("File not present in given location")
    try: 
        return at.Table.read(location)
    except: 
        raise NameError("File cannot be opened as an astropy table")

def _read_tf(trans_func_loc):
    """ Returns a 1D interpolated function of input """
    trans_func = _initialise(trans_func_loc)
    x_axis = np.array([i for i in trans_func['ws']])
    y_axis = np.array([j for j in trans_func['trans']])
    return interpolate.interp1d(x_axis, y_axis)


def _limits(band):
    """ Returns initialisation wavelength and termination wavelength """
    limits = {'g':(3990.0, 5440.0), 'r':(5440.0, 6960.0), 'r':(5440.0, 6960.0),
              'i':(6960.0, 8500.0), 'z':(8500.0, 9320.0), 'y':(9320.0, 10760.0)}
    return limits[band]
    
def _get_params(parent_dir, Lambda, _flux, band):
    """ Read and return req. parameters for convolution """
    init, ter = _limits(band)
    # _lambda contains Lambda values for particular band
    _lambda = [l for l in Lambda if l >= init and l<=ter]
    min_lambda = [i for i in _lambda if i>=init]
    init_flux = Lambda.index(min_lambda[0])
    _lambda = np.array(_lambda)
    # Interpolated T_band is present in T_band
    foo = _read_tf(parent_dir+band+".csv")
    T_band = np.array(foo(_lambda))
    flux = _flux[init_flux:init_flux+len(_lambda)]
    F_lambda = np.array([i for i in flux])
    
    return _lambda, T_band, F_lambda

def _integrate(function, Lambda):
    """ Performs integration by trapezoidal rule """
    # Variable definition for Area
    Area = 0

    # loop to calculate the area of each trapezoid and sum.
    for i in range(1, len(Lambda)):
        #the x locations of the left and right side of each trapezpoid
        x0 = Lambda[i]
        x1 = Lambda[i-1]

        """ 'dx' will be the variable width of the trapezoid """
        dx = abs(x1 - x0)
        
        """
        TRAPEZOIDAL RULE
        """
        """ Area of each trapezoid """
        Ai = dx * (function[i] + function[i-1])/ 2.
        """ Cumulative sum of the areas """
        Area += Ai
    
    return Area

def _calc_AB_magnitude(integral_1, integral_2):
    """ Calculating AB Magnitude from Integrals """
    temp = ((integral_1/integral_2)/(3 * 10**8))*(1e-17)*(10**-10)
    return -2.5*np.log10(temp) - 48.600

def _extrapolate(signal):
    """ Perform extrapolation on signal at first before anything else """
    _flux = signal['flux']
    _lambda = [l for l in 10**signal['loglam']]
    _flux = [f for f in _flux]
    
    if max(_lambda)<10760.0:
        step = abs(_lambda[-1] - _lambda[-2])
        add_lambda = range(max(_lambda), 10760, step)
        _lambda.extend(add_lambda)
        _flux.extend([_flux[-1]]*len(add_lambda))
    return _lambda, _flux
        
def get_AB_magnitudes(spectrum_loc):
    """ call get_params and get required values for that signal and all bands """
    parent_dir = ""
    AB_mag = []
    spectrum = _initialise(spectrum_loc)
    bands = ['g', 'r', 'i', 'z', 'y']
    # Extrapolate given sepctrum if required
    _lambda, _flux = _extrapolate(spectrum)
    for band in bands:
        Lambda, T_band, F_lambda = _get_params(parent_dir, _lambda, _flux,  band)
        """
        Numerator and Denominator Function
        """
        func1 = (Lambda * F_lambda * T_band)
        func2 = (T_band / Lambda) 
        
        Integral_1 = _integrate(func1, Lambda)
        Integral_2 = _integrate(func2, Lambda)
        
        if Integral_1 or Integral_2:
            AB_mag.append(_calc_AB_magnitude(Integral_1, Integral_2))
        else:
            """ Extrapolate and calculating AB Magnitude """
            AB_mag.append(0)
        
    specMag_g, specMag_r, specMag_i, specMag_z, specMag_y = AB_mag
    
    # Returns in order of bands provided
    return specMag_g, specMag_r, specMag_i, specMag_z, specMag_y