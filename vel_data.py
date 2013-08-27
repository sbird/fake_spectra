# -*- coding: utf-8 -*-
"""Plot the metal line data from Prochaska"""

import numpy as np
import matplotlib.pyplot as plt
import leastsq as ls

def load_data(zrange=None):
    """Load the data tables and concatenate them"""
    #data = np.loadtxt("vel_width_mtlkin_704_diff.dat")
    #vel_data = data[:,2]
    data2 = np.loadtxt("apj469315t2_mrt_mod.txt")
    redshift = data2[:,0] #np.concatenate([data[:,0], data2[:,0]])
    met = data2[:,3] #np.concatenate([data[:,1],data2[:,3]])
    vel_data = data2[:,5] #np.concatenate([data[:,2], data2[:,5]])
    if zrange != None:
        ind = np.where((redshift < zrange[0])*(redshift > zrange[1]))
        return (redshift[ind], met[ind], vel_data[ind])
    return (redshift, met, vel_data)

def pdf_with_error(vel_data, nv_table):
    """Plot a pdf of vel_data in nv_table bins, with error bars"""
    v_table=np.logspace(np.log10(np.min(vel_data)),np.log10(np.max(vel_data)),nv_table)
    center = np.array([(v_table[i]+v_table[i+1])/2. for i in range(0,np.size(v_table)-1)])

    nn = np.histogram(vel_data,v_table)[0]
    #Normalise so that integral in log space is unity.
    #This emulates np.histogram(np.log10(met), np.log10(bin),density=True)
    width = np.array([(-np.log10(v_table[i])+np.log10(v_table[i+1])) for i in range(0,np.size(v_table)-1)])

    #This is the avg. fraction of DLAs per unit absorption distance.
    #It is needed to match the normalisation of Pontzen 2008.
    #DLAdX = 0.065

    norm = width*np.size(vel_data)

    vels = nn / norm
    #Use poisson errors
    verr = (np.sqrt(nn))/norm
    #These v_table will have error bars that go to 0
    ind = np.where(nn == 1)
    verr[ind] *= (0.98)
    plt.errorbar(center,vels,xerr=[center-v_table[:-1],v_table[1:]-center],yerr=verr,fmt='.', color="black")

    plt.semilogx(center, vels,'o')
    return (center, vels,verr)

def plot_prochaska_2008_data(zrange = None, nv_table=7):
    """Plot a velocity width histogram from Prochaska 2008/2013"""
    (_, _, vel_data) = load_data(zrange)
    (center, vels, verr) = pdf_with_error(vel_data, nv_table)
    plt.xlim(10, 1000)
    return (center, vels,verr)

def plot_alpha_metal_data(zrange=None,nv_table=7):
    """
       Plot the metallicities from alpha peak mainly (Si and S) elements
       from Neeleman 2013
       redshift  = (4,3) will show only quasars between z=4 and z=3
    """
    (_, met, _) = load_data(zrange)
    return pdf_with_error(10**met, nv_table)

def plot_prochaska_2008_correlation(zrange = None, color="black"):
    """Plot the observed correlation between velocity widths and metallicity from Prochaska 2008"""
    (_, met, vel_data) = load_data(zrange)
    plt.loglog(vel_data, 10**met,'o',color=color)
    vel = np.log10(vel_data)
    (intercept, slope, var) = ls.leastsq(met,vel)
    print "obs fit: ",intercept, slope, np.sqrt(var)
    print "obs correlation: ",ls.pearson(met, vel,intercept, slope)
    print "obs kstest: ",ls.kstest(met, vel,intercept, slope)
    xx = np.logspace(np.min(met), np.max(met),15)
    plt.loglog(10**intercept*xx**slope, xx, color=color)

def plot_extra_stat_hist(stat=False,zrange=None, nv_table=11):
    """Plot a histogram of the mean-median statistic"""
    data2 = np.loadtxt("apj469315t2_mrt_mod.txt")
    redshift = data2[:,0] #np.concatenate([data[:,0], data2[:,0]])
    fmm = data2[:,6+int(stat)] #np.concatenate([data[:,1],data2[:,3]])
    if zrange != None:
        ind = np.where((redshift < zrange[0])*(redshift > zrange[1]))
        fmm = fmm[ind]

    v_table=np.linspace(0,1,nv_table)
    center = np.array([(v_table[i]+v_table[i+1])/2. for i in range(0,np.size(v_table)-1)])

    nn = np.histogram(fmm,v_table)[0]
    #Normalise so that integral in log space is unity.
    #This emulates np.histogram(np.log10(met), np.log10(bin),density=True)
    width = np.array([(-v_table[i]+v_table[i+1]) for i in range(0,np.size(v_table)-1)])

    norm = width*np.size(fmm)

    vels = nn / norm
    #Use poisson errors
    verr = (np.sqrt(nn))/norm
    plt.errorbar(center,vels,xerr=[center-v_table[:-1],v_table[1:]-center],yerr=verr,fmt='.', color="black")
    plt.plot(center, vels, 'o')
