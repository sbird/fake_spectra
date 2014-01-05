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

outdir = path.join(myname.base, "plots/checks/")
print "Plots at: ",outdir

def plot_vel_width_metcol(sim, snap):
    """Load a simulation and plot its velocity width"""
    halo = myname.get_name(sim)
    #Load from a save file only
    hspec = ps.PlottingSpectra(snap, halo)
    hspec.plot_vel_width("Si", 2, met_cut = None)
    hspec.plot_vel_width("Si", 2, met_cut = 1e13, color="blue")
    vel_data.plot_prochaska_2008_data()
    save_figure(path.join(outdir,"cosmo"+str(sim)+"_low_metals_z"+str(snap)))
    plt.clf()
    (vbin, vels1) = hspec.vel_width_hist("Si", 2, met_cut = None)
    (vbin, vels2) = hspec.vel_width_hist("Si", 2, met_cut = 1e13)
    mm = np.min((np.size(vels2), np.size(vels1)))
    plt.semilogx(vbin[:mm], vels2[:mm]/vels1[:mm], color="black")
    save_figure(path.join(outdir,"cosmo"+str(sim)+"_low_metals_rel_z"+str(snap)))
    plt.clf()

def plot_vel_width_SiII(sim, snap, ff=False):
    """
       Plot the change in velocity widths between the full calculation and
       setting n(Si+)/n(Si) = n(HI)/n(H)
    """
    halo = myname.get_name(sim, ff)
    #Load from a save file only
    hspec = ps.PlotHaloSpectra(snap, halo)
    hspec.plot_vel_width("Si", 2)
    hspecSi = ps.PlotHaloSpectra(snap, halo,savefile="SiHI_spectra.hdf5")
    hspecSi.plot_vel_width("Si", 2, color="blue")
    vel_data.plot_prochaska_2008_data()
    save_figure(path.join(outdir,"cosmo"+str(sim)+"_SiHI_z"+str(snap)))
    plt.clf()
    (vbin, vels1) = hspec.vel_width_hist("Si", 2)
    (vbin, vels2) = hspecSi.vel_width_hist("Si", 2)
    mm = np.min((np.size(vels2), np.size(vels1)))
    plt.semilogx(vbin[:mm], vels2[:mm]/vels1[:mm], color="black")
    save_figure(path.join(outdir,"cosmo"+str(sim)+"_SiHI_rel_z"+str(snap)))
    plt.clf()

def test_spec_resolution():
    """Plot the velocity widths for different spectral resolutions"""
    #Do spectral resolution test
    halo = myname.get_name(7)
    #Higher resolution spectrum
    hspec = ps.PlottingSpectra(3, halo, savefile="grid_spectra_DLA.hdf5")
    hspec2 = ps.PlottingSpectra(3, halo, None, None, savefile="grid_spectra_DLA_res.hdf5")
    plot_check(hspec,hspec2,"specres")

def test_pecvel():
    """Plot the velocity widths with and without peculiar velocities"""
    #Do spectral resolution test
    halo = myname.get_name(7)
    #Higher resolution spectrum
    hspec = ps.PlottingSpectra(3, halo, savefile="grid_spectra_DLA.hdf5")
    hspec2 = ps.PlottingSpectra(3, halo, None, None, savefile="grid_spectra_DLA_pecvel.hdf5")
    plot_check(hspec,hspec2,"pecvel")

def test_box_resolution():
    """Plot the velocity widths for different size boxes"""
    #Do spectral resolution test
    halo = myname.get_name(5)
    halo10 = myname.get_name(5,box=10)
    hspec = ps.PlottingSpectra(3, halo)
    hspec2 = ps.PlottingSpectra(3, halo10)
    plot_check(hspec,hspec2,"box")

    #Higher resolution spectrum
def plot_check(hspec, hspec2, ofile):
    """Plot two halos both absolutely and relatively"""
    hspec.plot_vel_width("Si",2, color="red")
    hspec2.plot_vel_width("Si", 2, color="blue")
    vel_data.plot_prochaska_2008_data()
    plt.xlim(1, 1000)
    save_figure(path.join(outdir,"cosmo_vel_width_"+ofile))
    plt.clf()
    (vbin,one) = hspec.vel_width_hist("Si",2)
    (vbin,two) = hspec2.vel_width_hist("Si",2)
    if np.size(one) != np.size(two):
        maxx = np.min([np.size(one),np.size(two)])
        plt.semilogx(vbin[:maxx],one[:maxx]/two[:maxx])
    else:
        plt.semilogx(vbin,one/two)
    save_figure(path.join(outdir,"cosmo_rel_vel_width_"+ofile))
    plt.clf()

if __name__ == "__main__":
    test_spec_resolution()
    test_box_resolution()
    test_pecvel()

#     plot_vel_width_SiII(0, 3)

    plot_vel_width_metcol(0,3)
