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
from save_figure import save_figure

base="/home/spb/scratch/Cosmo/"
outdir = base + "plots/"
print "Plots at: ",outdir

def plot_vel_width_sim(sim, snap, color="red", ff=False, HI_cut = None):
    """Load a simulation and plot its velocity width"""
    halo = "Cosmo"+str(sim)+"_V6"
    if ff:
        halo+="_512"
    #Load from a save file only
    hspec = ps.PlottingSpectra(snap, base+halo, None, None)
    hspec.plot_vel_width("Si", 2, color=color, HI_cut = HI_cut)

def plot_sep_frac(sim, snap):
    """Plot fraction of lines from separated halos"""
    halo = "Cosmo"+str(sim)+"_V6"
    #Load from a save file only
    hspec = ps.PlottingSpectra(snap, base+halo, None, None)
    hspec.plot_sep_frac()

def plot_rel_vel_width(sim1, sim2, snap, color="black"):
    """Load and make a plot of the difference between two simulations"""
    halo1 = "Cosmo"+str(sim1)+"_V6"
    halo2 = "Cosmo"+str(sim2)+"_V6"
    hspec1 = ps.PlotHaloSpectra(snap, base+halo1)
    (vbin, vels1) = hspec1.vel_width_hist("Si", 2)
    hspec1 = ps.PlotHaloSpectra(snap, base+halo2)
    (vbin, vels2) = hspec1.vel_width_hist("Si", 2)
    mm = np.min((np.size(vels2), np.size(vels1)))
    plt.semilogx(vbin[:mm], vels2[:mm]/vels1[:mm], color=color)

def plot_spectrum(sim, snap, num):
    """Plot a spectrum"""
    halo = "Cosmo"+str(sim)+"_V6"
    #Load from a save file only
    hspec = ps.PlottingSpectra(snap, base+halo, None, None)
    tau = hspec.get_tau("Si", 2)
    hspec.plot_spectrum(tau, num)

def plot_spectrum_density_velocity(sim, snap, num):
    """Plot a spectrum"""
    halo = "Cosmo"+str(sim)+"_V6"
    #Load from a save file only
    hspec = ps.PlotHaloSpectra(snap, base+halo)
    hspec.plot_spectrum_density_velocity("Si",2, num)

def plot_metallicity(sim, snap):
    """Plot a spectrum"""
    halo = "Cosmo"+str(sim)+"_V6"
    #Load from a save file only
    hspec = ps.PlottingSpectra(snap, base+halo, None, None, savefile="halo_spectra_DLA.hdf5")
    hspec.plot_metallicity()
    vel_data.plot_alpha_metal_data(nbins=12)
    save_figure(path.join(outdir,"cosmo"+str(sim)+"_metallicity_z"+str(snap)))
    plt.clf()
    hspec.plot_Z_vs_vel_width()
    vel_data.plot_prochaska_2008_correlation()
    save_figure(path.join(outdir,"cosmo"+str(sim)+"_correlation_z"+str(snap)))
    plt.clf()

if __name__ == "__main__":
    colors=["blue", "purple", "orange", "red"]
    reds = {54:4, 60:3, 68:2}
    plot_spectrum(0,60, 102)
    plt.xlim(1100, 1600)
    save_figure(path.join(outdir,"cosmo0_Si_spectrum"))
    plt.clf()

    plot_metallicity(0, 60)

    for ss in (0,1,2,3):
        plot_spectrum_density_velocity(ss,60, 25)
        save_figure(path.join(outdir,"cosmo"+str(ss)+"_Si_spectrum"))
        plt.clf()

    for ss in (0,1,2,3):
        plot_sep_frac(ss,60)
    save_figure(path.join(outdir,"cosmo_sep_frac_z3"))
    plt.clf()



    #Best-fit base model
    plot_vel_width_sim(0, 60, "blue", HI_cut = 10**20.3)
    vel_data.plot_prochaska_2008_data()
    plt.ylim(1e-5, 2e-2)
    plt.xlim(10, 1000)
    save_figure(path.join(outdir,"cosmo_vel_width_z3"))
    plt.clf()

    for sp in (60,68):
        #The vel widths for different simulations
        for ss in (3,0):
            plot_vel_width_sim(ss, sp, colors[ss], HI_cut = 10**17)

        vel_data.plot_prochaska_2008_data()
        plt.ylim(1e-5, 2e-2)
        plt.xlim(10, 1000)
        save_figure(path.join(outdir,"cosmo_feedback_z"+str(reds[sp])))
        plt.clf()

    #A plot of the redshift evolution
    zz = [54,60,68]
    for ii in (0,1,2):
        plot_vel_width_sim(0, zz[ii], color=colors[ii])

    vel_data.plot_prochaska_2008_data()
    plt.title("Velocity Widths at z=4-2")
    save_figure(path.join(outdir,"cosmo_vel_width_zz"))
    plt.clf()

    plot_rel_vel_width(1,0,60)
    plot_rel_vel_width(2,0,60, color="grey")
    #Something odd about this one...
    plot_rel_vel_width(3,0,60, color="red")

    save_figure(path.join(outdir,"cosmo_rel_vel_width_z3"))
    plt.clf()
