#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
Project Details:  
    outlier_detect.py 
    Detects and removes outliers using diffrent methods
    Choice of algorithm is upto the user

Created on Thu Aug 24 07:49:57 2017

__author__      = "nnarenraju"
__copyright__   = "Copyright 2017, AGN_classification"
__credits__     = "nnarenraju"
__license__     = "Apache License 2.0"
__version__     = "1.0.1"
__maintainer__  = "nnarenraju"
__email__       = "nnarenraju@gmail.com"
__status__      = "inUsage"

Github Repository: "https://github.com/nnarenraju/ASIAA_SSP"

"""

import os
import csv
import numpy as np
from scipy import stats
import astropy.table as at
import matplotlib.font_manager
import matplotlib.pyplot as plt

from sklearn import svm
from sklearn.covariance import EllipticEnvelope
from sklearn.ensemble import IsolationForest

# Initialise random
rng = np.random.RandomState(42)

def _initialise(location):
    """ Returns an astropy table if possible """
    if not os.path.exists(location):
        raise NameError("File not present in given location")
    try: 
        return at.Table.read(location)
    except: 
        raise NameError("File cannot be opened as an astropy table")
        

def _classifiers(n_samples):
    """ Returns connfiguration of classifiers used """
    classifiers = {
            "One-Class SVM": svm.OneClassSVM(kernel="rbf", gamma=0.1),
            "Robust covariance": EllipticEnvelope(),
            "Isolation Forest": IsolationForest(max_samples=n_samples, random_state=rng)}
    
    return classifiers


def _plot_contour(X, clf, n_errors, clf_name, outliers_fraction, i):
    """ Plots contours with given dataset and detected outliers """
    # Sanity Check 0
    if len(X[0])!=2:
        raise ValueError("Input data dimension incorrect, cannot be plotted on 2D graph")
    
    xx, yy = np.meshgrid(np.linspace(-7, 7, 100), np.linspace(-7, 7, 100))
    scores_pred = clf.decision_function(X)
    threshold = stats.scoreatpercentile(scores_pred, 100 * outliers_fraction)
    
    # plot the levels lines and the points
    Z = clf.decision_function(np.c_[xx.ravel(), yy.ravel()])
    Z = Z.reshape(xx.shape)
        
    subplot = plt.subplot(2, 2, i + 1)
    subplot.contourf(xx, yy, Z, levels=np.linspace(Z.min(), threshold, 7), cmap=plt.cm.Blues_r)
    a = subplot.contour(xx, yy, Z, levels=[threshold], linewidths=2, colors='red')
        
    subplot.contourf(xx, yy, Z, levels=[threshold, Z.max()], colors='orange')
    b = subplot.scatter(X[:-n_errors, 0], X[:-n_errors, 1], c='white', s=20, edgecolor='k')
    c = subplot.scatter(X[-n_errors:, 0], X[-n_errors:, 1], c='black', s=20, edgecolor='k')
        
    subplot.axis('tight')
    subplot.legend(
        [a.collections[0], b, c],
        ['learned decision function', 'true inliers', 'true outliers'],
        prop=matplotlib.font_manager.FontProperties(size=10),
        loc='lower right')
    subplot.set_xlabel("%d. %s (errors: %d)" % (i + 1, clf_name, n_errors))
    subplot.set_xlim((-7, 7))
    subplot.set_ylim((-7, 7))
    plt.subplots_adjust(0.04, 0.1, 0.96, 0.94, 0.1, 0.26)
    plt.suptitle("Outlier detection")


# Fit the problem with varying cluster separation
def _ML_detect(data, plot_contour=False):
    """ Detects outliers using three different algorithms """
    np.random.seed(42)
    # Input data
    X, n_samples = data, len(data)
    
    # Classifiers
    classifiers = _classifiers(n_samples)
    
    def _get_dupe(seq, item):
        """ Returns locations of item in given list if duplicated/non-duplicated """
        start_at = -1
        locs = []
        while True:
            try:
                loc = seq.index(item, start_at+1)
            except ValueError:
                break
            else:
                locs.append(loc)
                start_at = loc
        return locs

    # Fit the model
    plt.figure(figsize=(9, 7))
    for i, (clf_name, clf) in enumerate(classifiers.items()):
        # fit the data and tag outliers
        clf.fit(X)
        y_pred = clf.predict(X)
        outliers = _get_dupe(y_pred, item=-1)
        
        """ Find the outlier number and outlier fraction to plot contour and have cool title """
        n_errors = len([i for i in y_pred if i==-1])
        outliers_fraction = float(n_errors)/float(n_samples)
        
        if plot_contour:
            _plot_contour(X, clf, n_errors, clf_name, outliers_fraction, i)

        return outliers


def _get_outliers(table, zrange, method="MAD", print_details=False):
    """ Get the specMag data and segregate accordingly """
    
    # Sanity Check 0
    init, ter = zrange
    if init<=0.0 or ter>=1.0:
        # Checks whether redshift is within known limits
        raise ValueError("Discrepancy in redshift range")
    
    # Getting redshift within specified limits
    _z = table[table['z'] > init] 
    table = _z[_z['z'] < ter]
     
    """       
    --> Abnormally noisy data will have higher magnitudes compared to other spectrums.
    
    --> As the data distribution is highly variant between different ranges of redshift,
        with few ranges having very few data instances, will be profoundly affected by the
        prescence of outliers and thus needs to be removed.
        
    --> Taking the sum of all convolved magnitudes and distinguishing data wrt them, in
        different ranges of redshift, will remove potential features that would prove 
        detrimental to the regression process.
    """
    
    ML_data = []
    for i in table:
        ML_data.append(i['specMag_g']+i['specMag_r']+i['specMag_i']+i['specMag_z']+i['specMag_y'])
    MAD_data = np.array(ML_data)

    if method=="MAD":
        outliers = _MAD_detect(MAD_data)
    elif method=="ML":
        outliers = _ML_detect(ML_data)
    
    if print_details:
        print "Z-Range: {0} to {1}".format(zrange[0], zrange[1])
        print "Length of z_range instances: ", len(MAD_data)
        print "Number of outliers in given z-range: ", len(outliers)
        print "Fraction of outliers present in given z-range: ", len(outliers)/len(MAD_data)
    
    return outliers


def _remove_outliers(outliers):
    """ Removes outliers and creates outlier free data as a CSV """
    outliers = list(outliers)    
    temp = table[:]
    indices = []
    for i in temp:
        target = i['specMag_g']+i['specMag_r']+i['specMag_i']+\
                 i['specMag_z']+i['specMag_y']
        if  target in outliers:
            indices.append(outliers.index(target))
    temp.remove_rows(indices)
    
    outlier_free_table = temp
    with open('outlier_free_table.csv', 'wb') as f:
        writer = csv.writer(f)
        writer.writerows(outlier_free_table)
    
    print "Outliers have been succesfully removed and a refined table has been saved as a CSV"


def _MAD_detect(x, printing=False):
    outliers = x[_mad_based_outlier(x)]
    if printing:
        print outliers
    return list(outliers)


def _mad_based_outlier(points, thresh=20):
    if len(points.shape) == 1:
        points = points[:,None]
    median = np.median(points, axis=0)
    diff = np.sum((points - median)**2, axis=-1)
    diff = np.sqrt(diff)
    med_abs_deviation = np.median(diff)
    modified_z_score = 0.6745 * diff / med_abs_deviation
    return modified_z_score > thresh


def _plot(x, y):
    """ Plots the distribution of data without outliers """
    #### Plotting Data Distribution
    plt.figure(figsize=(16,9))
    plt.grid(True)
    plt.title("Distribution of Data in Input Dataset", fontsize=18)
    plt.xticks(fontsize=14)  
    plt.yticks(fontsize=12)  
    plt.xlabel("Redshift (z) = [x_value*0.01] + 0.02", fontsize=15)
    plt.ylabel("Number of occurences", fontsize=15)  
    plt.bar(x, y, align='center', color="#3F5D7D")
    plt.savefig('Data_distribution.png')
    plt.show()
    

if __name__ == "__main__":
    
    #Initial range
    z_ranges = [(0.20, 0.29999), (0.30, 0.39999), (0.40, 0.49999), (0.50, 0.59999),
               (0.60, 0.69999), (0.70, 0.79999), (0.80, 0.89999), (0.90, 0.99999) ]
    
    # Initialse the input data
    location = 'C:\Users\Admin\Desktop\dataset_ML.csv'
    table = _initialise(location)
    
    outliers = []
    for z_range in z_ranges:
        outliers.append(_get_outliers(table, z_range, method="MAD", print_details=True))

    print "Total Number of outliers in given dataset: ", len(outliers)
    print "Fraction of outliers present: ", len(outliers)/len(table)
    print
    
    try:
        _plot(x, y)
    except:
        print "Plotting function under construction"
