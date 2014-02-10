#!/usr/bin env python
# -*- coding: utf-8 -*-
"""Make some plots of the velocity widths from the cosmo runs"""

import matplotlib
matplotlib.use('PDF')

import matplotlib.pyplot as plt

import plot_spectra as ps
import vel_data
import os.path as path
import numpy as np
import myname
from save_figure import save_figure

outdir = path.join(myname.base, "plots/")
print "Plots at: ",outdir
zrange = {1:(7,3.5), 3:(7,0), 5:(2.5,0)}
#Colors and linestyles for the simulations
colors = {0:"red", 1:"purple", 2:"cyan", 3:"green", 4:"gold", 5:"orange", 7:"blue", 6:"grey"}
lss = {0:"--",1:":", 2:":",3:"-.", 4:"--", 5:"-",6:"--",7:"-"}
labels = {0:"ILLUS",1:"HVEL", 2:"HVNOAGN",3:"NOSN", 4:"WMNOAGN", 5:"MVEL",6:"METAL",7:"2xUV"}

hspec_cache = {}

def get_hspec(sim, snap):
    """Get a spectra object, possibly from the cache"""
    halo = myname.get_name(sim, True)
    #Load from a save file only
    try:
        hspec = hspec_cache[(sim, snap)]
    except KeyError:
        hspec = ps.PlottingSpectra(snap, halo, label=labels[sim])
        hspec_cache[(sim, snap)] = hspec
    return hspec

def plot_vel_width_sim(sim, snap, color="red", HI_cut = None):
    """Load a simulation and plot its velocity width"""
    hspec = get_hspec(sim, snap)
    hspec.plot_vel_width("Si", 2, color=color, HI_cut = HI_cut)

def plot_sep_frac(sim, snap):
    """Plot fraction of lines from separated halos"""
    hspec = get_hspec(sim, snap)
    hspec.plot_sep_frac()

def plot_spectrum(sim, snap, num, subdir=""):
    """Plot a spectrum"""
    hspec = get_hspec(sim, snap)
    tau = hspec.get_observer_tau("Si", 2)
    ind = hspec.get_filt("Si", 2)
    (low,high, offset) = hspec.find_absorber_width("Si",2, minwidth=0)
    tau_l = np.roll(tau[ind][num], offset[ind][num])
    hspec.plot_spectrum(tau_l[low[ind][num]:high[ind][num]])
    save_figure(path.join(outdir,"spectra/"+subdir+"cosmo"+str(sim)+"_Si_"+str(num)+"_spectrum"))
    plt.clf()
    tau = hspec.get_tau("Si", 2,4)
    tau_l = np.roll(tau[ind][num], offset[ind][num])
    hspec.plot_spectrum(tau_l[low[ind][num]:high[ind][num]])
    save_figure(path.join(outdir,"spectra/"+subdir+"cosmo"+str(sim)+"_Si_"+str(num)+"_1260_spectrum"))
    plt.clf()

def plot_spectrum_max(sim, snap):
    """Plot spectrum with max vel width"""
    hspec = get_hspec(sim, snap)
    vels = hspec.vel_width("Si",2)
    ind = hspec.get_filt("Si", 2)
    num = np.where(vels[ind] == np.max(vels[ind]))[0][0]
    plot_spectrum(sim, snap, num, "max/")
    num = np.where(vels[ind] > 500)[0]
    for nn in num:
        plot_spectrum(sim, snap, nn, "max/")


def plot_spectrum_density_velocity(sim, snap, num):
    """Plot a spectrum"""
    hspec = get_hspec(sim, snap)
    hspec.plot_spectrum_density_velocity("Si",2, num)
    save_figure(path.join(outdir,"spectra/cosmo"+str(sim)+"_tdv_Si_spectrum"))
    plt.clf()

def plot_metallicity(sims, snap):
    """Plot metallicity, vel width, their correlation and the extra statistics"""
    out = "cosmo_metallicity_z"+str(snap)
    for sim in sims:
        hspec = get_hspec(sim, snap)
        hspec.plot_metallicity(color=colors[sim], ls=lss[sim])
    vel_data.plot_alpha_metal_data(zrange[snap])
    plt.legend(loc=2,ncol=3)
    plt.ylim(0,2)
    save_figure(path.join(outdir,out))
    plt.clf()

def plot_met_corr(sims,snap):
    """Plot metallicity velwidth correlations"""
    for sim in sims:
        print "Met vel corr. Simulation ",sim
        out = "cosmo"+str(sim)+"_correlation_z"+str(snap)
        hspec = get_hspec(sim, snap)
        hspec.plot_Z_vs_vel_width(color=colors[sim])
        vel_data.plot_prochaska_2008_correlation(zrange[snap])
        save_figure(path.join(outdir,out))
        plt.clf()
        (_, met, vels) = vel_data.load_data()
        print "KS test is : ",hspec.kstest(10**met, vels)

def plot_vel_width_sims(sims, snap):
    """Plot velocity widths for a series of simulations"""
    vel_data.plot_prochaska_2008_data(zrange[snap])
    for sss in sims:
        #Make abs. plot
        hspec = get_hspec(sss, snap)
        hspec.plot_vel_width("Si", 2, color=colors[sss], ls=lss[sss])
    plt.ylim(0,2)
    plt.xlim(10,1000)
    plt.legend(loc=2,ncol=3)
    save_figure(path.join(outdir,"cosmo_vel_width_z"+str(snap)))
    plt.clf()

def plot_rel_vel_width(sims, snap):
    """Plot velocity widths relative to simulation 7"""
    hspec = get_hspec(7, snap)
    (vbin, vels7) = hspec.vel_width_hist("Si", 2)
    #Make rel plot
    for sss in sims:
        hspec = get_hspec(sss, snap)
        (vbin, vel) = hspec.vel_width_hist("Si", 2)
        mm = np.min([np.size(vel), np.size(vels7)])
        plt.semilogx(vbin[:mm], vel[:mm]/vels7[:mm], color=colors[sss],ls=lss[sss])
    plt.xlim(10, 1000)
    save_figure(path.join(outdir,"cosmo_rel_vel_z"+str(snap)))
    plt.clf()

def plot_mean_median(sims, snap):
    """Plot mean-median statistic for all sims on one plot"""
    #Plot extra statistics
    for sss in sims:
        hspec = get_hspec(sss, snap)
        hspec.plot_f_meanmedian("Si", 2, color=colors[sss], ls=lss[sss])
    vel_data.plot_extra_stat_hist(False,zrange[snap])
    plt.ylim(0,3)
    plt.legend(loc=2,ncol=3)
    save_figure(path.join(outdir,"cosmo_mean_median_z"+str(snap)))
    plt.clf()

def plot_f_peak(sims, snap):
    """Plot peak statistic for all sims on one plot"""
    for sss in sims:
        hspec = get_hspec(sss, snap)
        hspec.plot_f_peak("Si", 2, color=colors[sss], ls=lss[sss])
    plt.legend(loc=2,ncol=3)
    vel_data.plot_extra_stat_hist(True,zrange[snap])
    plt.ylim(0,3)
    save_figure(path.join(outdir,"cosmo_peak_z"+str(snap)))
    plt.clf()

def plot_vel_widths_cloudy():
    """Plot some velocity width data for different cloudy models"""
    #Load sims
    hspec0 = ps.PlottingSpectra(3, myname.get_name(0, True))
    hspec1 = ps.PlottingSpectra(3, myname.get_name(0,True), savefile="rand_spectra_DLA_fancy_atten.hdf5")
    #Make abs. plot
    hspec0.plot_vel_width("Si", 2, color="blue", ls="--")
    hspec1.plot_vel_width("Si", 2, color="red", ls="-")
    vel_data.plot_prochaska_2008_data()
    save_figure(path.join(outdir,"cosmo_feedback_cloudy_z3"))
    plt.clf()
    #Make rel plot
    (vbin, vels0) = hspec0.vel_width_hist("Si", 2)
    (vbin, vels2) = hspec1.vel_width_hist("Si", 2)
    mm = np.min((np.size(vels2),np.size(vels0)))
    plt.semilogx(vbin[:mm], vels0[:mm]/vels2[:mm], color="blue",ls="-")
    plt.xlim(1, 1000)
    save_figure(path.join(outdir,"cosmo_rel_vel_cloudy_z3"))
    plt.clf()

def plot_vel_redshift_evo(sim):
    """Plot the evolution with redshift of a simulation"""
    halo = myname.get_name(sim, True)
    vels = {}
    for snap in (1,3,5):
        hspec0 = ps.PlottingSpectra(snap, halo)
        (vbin, vels[snap]) = hspec0.vel_width_hist("Si", 2)
    mm = np.min([np.size(vel) for vel in vels.values()])
    #Normalised by z=3
    plt.semilogx(vbin[:mm], vels[3][:mm]/vels[5][:mm], color="black",ls="--")
    plt.semilogx(vbin[:mm], vels[1][:mm]/vels[5][:mm], color="grey",ls="-")
    plt.xlim(1, 1000)
    plt.ylim(0,2)
    save_figure(path.join(outdir,"cosmo_"+str(sim)+"_zz_evol"))
    plt.clf()

if __name__ == "__main__":
#     plot_vel_widths_cloudy()

    simlist = (0,1,3,7) #range(8)
    for ss in simlist:
        plot_spectrum_density_velocity(ss,3, 15)
        for nn in (272,350,457,1030,1496,2030,3333):
            plot_spectrum(ss,3, nn)
        plot_spectrum_max(ss,3)

    for zz in (1, 3, 5):
        plot_met_corr(simlist,zz)
        plot_metallicity(simlist, zz)
        plot_vel_width_sims(simlist, zz)
        plot_mean_median(simlist, zz)
        plot_f_peak(simlist, zz)

    for ss in (0,1,3,7):
        plot_sep_frac(ss,3)
    save_figure(path.join(outdir,"cosmo_sep_frac_z3"))
    plt.clf()

#     for ss in (0,1,2,3,4):
#         plot_vel_redshift_evo(ss)
