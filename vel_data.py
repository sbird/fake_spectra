# -*- coding: utf-8 -*-
"""Plot the metal line data from Prochaska"""

import numpy as np
import matplotlib.pyplot as plt
import leastsq as ls

def plot_prochaska_2008_data():
    """Plot a velocity width histogram from Prochaska 2008"""
    data = np.loadtxt("vel_width_mtlkin_704.dat")
    vel_data = data[:,2]

    #These are the same bins Tescari uses
    v_table = np.array([12,40,80,120,160,200,240,320,390,590,800])
    center = np.array([(v_table[i]+v_table[i+1])/2. for i in range(0,np.size(v_table)-1)])

    nn = np.histogram(vel_data,v_table)[0]
    #Bin width
    width = np.array([v_table[i+1]-v_table[i] for i in range(0,np.size(v_table)-1)])

    #This is the avg. fraction of DLAs per unit absorption distance.
    #It is needed to match the normalisation of Pontzen 2008.
    #DLAdX = 0.065

    norm = width* (1.*np.size(vel_data))

    vels = nn / norm
    #Use poisson errors
    verr = (np.sqrt(nn))/norm
    #These bins will have error bars that go to 0
    ind = np.where(nn == 1)
    verr[ind] *= (0.98)
    plt.errorbar(center,vels,xerr=[center-v_table[:-1],v_table[1:]-center],yerr=verr,fmt='.', color="black")

    plt.loglog(center, vels,'o')
    return (center, vels,verr)

def plot_alpha_metal_data(redshift=None,include_limits=False,color="black",nbins=9):
    """
       Plot the metallicities from alpha peak mainly (Si and S) elements
       from Prochaska 2007 (astro-ph/0702325).
       Some hints of bimodality here!
       Include also the metallicity histogram from Prochaska 2008 (astro-ph/0703701)
       redshift  = (4,3) will show only quasars between z=4 and z=3
    """

    #Columns: z_em      z_abs         NHI       error       f_a       α /H      error       f_zn     Zn/H       error       f_Fe    Fe/H       error
    #f key:  f_a: 0 = no measurement (set to -30). 1 = Si 2 = Si lower limit (error set to 0) 3 = Si upper limit 4 = S/H 5 = O/H 13 = Si + S limits
    # f_Zn:  0=No measurement; 1=Zn measurement; 2=Zn lower limit; 3=Zn upper limit
    # f_Fe: 0=No measurement; 1=Fe measurement; 2=Fe lower limit; 3=Fe upper limit; 4=[Ni/H]-0.1dex; 5=[Cr/H] - 0.2dex; 6=[Al/H]; 11-16=Fe, Ni, Cr, Al limits; 25=Mean of Fe limits.
    data = np.loadtxt("met_data_pro_2007.txt")
    #Select redshift
    if redshift != None:
        ind = np.where(np.logical_and(data[:,1] < redshift[0],data[:,1] > redshift[1]))
        data = data[ind]
    #Select on f_a: we only want firm measurements
    if not include_limits:
        ind = np.where(np.logical_or(data[:,4] == 1, data[:,4] > 3))
    else:
        ind = np.where(data[:,4] > 0)
    data2 = np.loadtxt("vel_width_mtlkin_704.dat")

    if redshift != None:
        zzind = np.where(np.logical_and(data2[:,0] < redshift[0],data2[:,0] > redshift[1]))
        data2 = data2[zzind]

    met = np.concatenate([data2[:,1], data[ind,5][0]])
    #met = 10**data[:,5]
    #Andrew's binning
    #nbins = 9
    big = np.max(np.concatenate([data2[:,1], data[ind,5][0]+data[ind,6][0]]))
    small = np.min(np.concatenate([data2[:,1], data[ind,5][0]-data[ind,6][0]]))
    bins=np.logspace(small,big,nbins)
    mbin = np.array([(bins[i]+bins[i+1])/2. for i in range(0,np.size(bins)-1)])
    nn = np.histogram(10**met,bins)[0]
    #Normalise so that integral in log space is unity.
    #This emulates np.histogram(np.log10(met), np.log10(bin),density=True)
    width = np.array([(-np.log10(bins[i])+np.log10(bins[i+1])) for i in range(0,np.size(bins)-1)])
    norm = width*np.size(met)
    hist = nn / norm
    merr = (np.sqrt(nn))/norm
    #These bins will have error bars that go to 0
    ind = np.where(nn == 1)
    merr[ind] *= (0.98)
    plt.errorbar(mbin,hist,xerr=[mbin-bins[:-1],bins[1:]-mbin],yerr=merr,fmt='.', color=color)
    plt.semilogx(mbin, hist,'o')

def plot_prochaska_2008_correlation(color="black"):
    """Plot the observed correlation between velocity widths and metallicity from Prochaska 2008"""
    data = np.loadtxt("vel_width_mtlkin_704.dat")
    met = data[:,1]
    vel = np.log10(data[:,2])
    plt.loglog(10**vel, 10**met,'o',color=color)
    (intercept, slope, var) = ls.leastsq(met,vel)
    print "obs fit: ",intercept, slope, np.sqrt(var)
    print "obs correlation: ",ls.pearson(met, vel,intercept, slope)
    xx = np.logspace(np.min(met), np.max(met),15)
    plt.loglog(10**intercept*xx**slope, xx, color=color)
