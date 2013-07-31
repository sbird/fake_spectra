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
zrange = {1:(7,3.5), 3:None, 5:(2.5,0)}
colors = {0:"red", 1:"purple", 2:"blue", 3:"green", 4:"orange"}
lss = {0:"--",1:":", 2:"-",3:"-.", 4:"-"}

def plot_vel_width_sim(sim, snap, color="red", ff=True, HI_cut = None):
    """Load a simulation and plot its velocity width"""
    halo = myname.get_name(sim, ff)
    #Load from a save file only
    hspec = ps.PlottingSpectra(snap, halo, None, None)
    hspec.plot_vel_width("Si", 2, color=color, HI_cut = HI_cut)

def plot_sep_frac(sim, snap):
    """Plot fraction of lines from separated halos"""
    halo = myname.get_name(sim, True)
    #Load from a save file only
    hspec = ps.PlottingSpectra(snap, halo, None, None)
    hspec.plot_sep_frac()

def plot_rel_vel_width(sim1, sim2, snap, color="black"):
    """Load and make a plot of the difference between two simulations"""
    halo1 = myname.get_name(sim1, True)
    halo2 = myname.get_name(sim2, True)
    hspec1 = ps.PlottingSpectra(snap, halo1)
    (vbin, vels1) = hspec1.vel_width_hist("Si", 2)
    hspec1 = ps.PlottingSpectra(snap, halo2)
    (vbin, vels2) = hspec1.vel_width_hist("Si", 2)
    mm = np.min((np.size(vels2), np.size(vels1)))
    plt.semilogx(vbin[:mm], vels2[:mm]/vels1[:mm], color=color)

def plot_spectrum(sim, snap, num):
    """Plot a spectrum"""
    halo = myname.get_name(sim, True)
    #Load from a save file only
    hspec = ps.PlottingSpectra(snap, halo, None, None)
    tau = hspec.get_observer_tau("Si", 2, num)
    hspec.plot_spectrum(tau)
    plt.xlim(1000,1500)

    save_figure(path.join(outdir,"cosmo"+str(sim)+"_Si_spectrum"))
    plt.clf()
    vels = hspec.vel_width(hspec.get_observer_tau("Si",2))
    ind = np.where(vels == np.max(vels[hspec.get_filt("Si",2)]))[0][0]
    tau2 = hspec.get_observer_tau("Si",2,ind)
    hspec.plot_spectrum(tau2)
    save_figure(path.join(outdir,"cosmo"+str(sim)+"_maxv_Si_spectrum"))
    plt.clf()


def plot_spectrum_density_velocity(sim, snap, num):
    """Plot a spectrum"""
    halo = myname.get_name(sim, True)
    #Load from a save file only
    hspec = ps.PlottingSpectra(snap, halo)
    hspec.plot_spectrum_density_velocity("Si",2, num)
    save_figure(path.join(outdir,"cosmo"+str(sim)+"_tdv_Si_spectrum"))
    plt.clf()

def plot_metallicity(sim, snap, ff=True):
    """Plot a spectrum"""
    halo = myname.get_name(sim, ff)
    out = "cosmo"+str(sim)+"_metallicity_z"+str(snap)
    if ff:
        out+="_512"
    #Load from a save file only
    hspec = ps.PlottingSpectra(snap, halo, None, None)
    hspec.plot_metallicity()
    vel_data.plot_alpha_metal_data(zrange[snap])
    save_figure(path.join(outdir,out))
    plt.clf()
    out = "cosmo"+str(sim)+"_correlation_z"+str(snap)
    if ff:
        out+="_512"
    hspec.plot_Z_vs_vel_width()
    vel_data.plot_prochaska_2008_correlation()
    save_figure(path.join(outdir,out))
    plt.clf()
    (redshift, met, vels) = vel_data.load_data()
    print "KS test is : ",hspec.kstest(10**met, vels)

def plot_vel_widths_sims(snap):
    """Plot some velocity width data at a particular redshift"""
    #Load sims
    hspec={}
    for sss in (0,1, 2, 3,4):
        hspec[sss] = ps.PlottingSpectra(snap, myname.get_name(sss,True))
        #Make abs. plot
        hspec[sss].plot_vel_width("Si", 2, color=colors[sss], ls=lss[sss])
    if snap == 3:
        vel_data.plot_prochaska_2008_data(zrange[snap], 11)
    else:
        vel_data.plot_prochaska_2008_data(zrange[snap])
    plt.ylim(0,1.6)
    save_figure(path.join(outdir,"cosmo_feedback_z"+str(snap)))
    plt.clf()
    vels={}
    #Make rel plot
    for sss in (0,2,3):
        (vbin, vels[sss]) = hspec[sss].vel_width_hist("Si", 2)
    mm = np.min([np.size(vel) for vel in vels.values()])
    plt.semilogx(vbin[:mm], vels[0][:mm]/vels[2][:mm], color=colors[0],ls=lss[0])
    plt.semilogx(vbin[:mm], vels[3][:mm]/vels[2][:mm], color=colors[3],ls=lss[3])
    plt.xlim(1, 1000)
    save_figure(path.join(outdir,"cosmo_rel_vel_z"+str(snap)))
    plt.clf()

def plot_vel_widths_res(snap):
    """Plot some velocity width data at a particular redshift"""
    #Load sims
    hspec0 = ps.PlottingSpectra(snap, myname.get_name(0, False))
    hspec512 = ps.PlottingSpectra(snap, myname.get_name(0,True))
    #Make abs. plot
    hspec0.plot_vel_width("Si", 2, color="blue", ls="--")
    hspec512.plot_vel_width("Si", 2, color="red", ls="-")
    vel_data.plot_prochaska_2008_data(zrange[snap])
    save_figure(path.join(outdir,"cosmo_feedback_res_z"+str(snap)))
    plt.clf()
    #Make rel plot
    (vbin, vels0) = hspec0.vel_width_hist("Si", 2)
    (vbin, vels2) = hspec512.vel_width_hist("Si", 2)
    mm = np.min((np.size(vels2),np.size(vels0)))
    plt.semilogx(vbin[:mm], vels0[:mm]/vels2[:mm], color="blue",ls="-")
    plt.xlim(1, 1000)
    save_figure(path.join(outdir,"cosmo_rel_vel_res_z"+str(snap)))
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
#     colors=["blue", "purple", "orange", "red"]

    for ss in (0,1,2,3,4):
        print "Metallicity Simulation",ss
        for zz in (1,3,5):
            plot_metallicity(ss, zz,True)

    for ss in (0,1,2,3,4):
        plot_spectrum_density_velocity(ss,3, 15)
        plot_spectrum(ss,3, 457)
    plot_spectrum(2,3, 272)

    for ss in (0,1,2,3,4):
        plot_sep_frac(ss,3)
    save_figure(path.join(outdir,"cosmo_sep_frac_z3"))
    plt.clf()

    for zz in (1, 3, 5):
        plot_vel_widths_sims(zz)
#         plot_vel_widths_res(zz)
#

    for ss in (0,1,2,3,4):
        plot_vel_redshift_evo(ss)
