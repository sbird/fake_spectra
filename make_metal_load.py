#!/usr/bin env python
# -*- coding: utf-8 -*-
"""Make some plots of the velocity widths from the cosmo runs"""

import matplotlib
matplotlib.use('PDF')

import matplotlib.pyplot as plt

import plot_spectra as ps
import vel_data
import os.path as path
from save_figure import save_figure

base="/home/spb/scratch/MetalLoad/"
outdir = base + "plots/"
print "Plots at: ",outdir
zrange = {17:(7,3.5), 25:None, 36:(2.5,0)}

def plot_metallicity(sim, snap):
    """Plot a spectrum"""
    halo = "modified_128_a"+str(sim)+"_b1"
    out = "cosmo"+str(sim)+"_metallicity_z"+str(snap)
    #Load from a save file only
    hspec = ps.PlottingSpectra(snap, base+halo, None, None)
    hspec.plot_metallicity()
    vel_data.plot_alpha_metal_data(zrange[snap])
    save_figure(path.join(outdir,out))
    plt.clf()
    out = "cosmo"+str(sim)+"_correlation_z"+str(snap)
    hspec.plot_Z_vs_vel_width()
    vel_data.plot_prochaska_2008_correlation()
    save_figure(path.join(outdir,out))
    hspec.save_file()
    plt.clf()

def plot_vel_widths_sims(sim):
    """Plot some velocity width data at a particular redshift"""
    #Load sims
    halo = "modified_128_a"+str(sim)+"_b1"
    hspec0 = ps.PlottingSpectra(17, base+halo)
    hspec2 = ps.PlottingSpectra(25, base+halo)
    hspec3 = ps.PlottingSpectra(36, base+halo)
    #Make abs. plot
    hspec0.plot_vel_width("Si", 2, color="blue", ls="-")
    hspec2.plot_vel_width("Si", 2, color="orange", ls="--")
    hspec3.plot_vel_width("Si", 2, color="red", ls="-.")
    vel_data.plot_prochaska_2008_data()
    save_figure(path.join(outdir,"cosmo_velw_"+str(sim)+"_metal_zz"))
    plt.clf()

def plot_vel_widths_metal():
    """Plot some velocity width data at a particular redshift"""
    #Load sims
    colorss={1:"blue", 4:"purple", 3:"orange", 5:"red"}
    for sim in (1,4,5):
        halo = "modified_128_a"+str(sim)+"_b1"
        hspec0 = ps.PlottingSpectra(25, base+halo)
        #Make abs. plot
        hspec0.plot_vel_width("Si", 2, color=colorss[sim], ls="-")
    vel_data.plot_prochaska_2008_data()
    save_figure(path.join(outdir,"cosmo_velw_metal_load"))
    plt.clf()

if __name__ == "__main__":
    colors=["blue", "purple", "orange", "red"]
    for zz in (17,25,36):
        plot_metallicity(1, zz)

    plot_vel_widths_sims(1)

    for ss in (4,5):
        plot_metallicity(ss, 25)

    plot_vel_widths_metal()


