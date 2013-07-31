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

def plot_vel_width_metcol(sim, snap, ff=False):
    """Load a simulation and plot its velocity width"""
    halo = myname.get_name(sim, ff)
    #Load from a save file only
    hspec = ps.PlottingSpectra(snap, halo, None, None)
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

def plot_vel_width_DLA(sim, snap, ff=False):
    """Plot the effect of the HI cut"""
    halo = myname.get_name(sim, ff)
    #Load from a save file only
    hspec = ps.PlottingSpectra(snap, halo, None, None)
    hspec.plot_vel_width("Si", 2, color="red", HI_cut = None)
    hspec.plot_vel_width("Si", 2, color="blue", HI_cut = 10**20.3)
    vel_data.plot_prochaska_2008_data()
    save_figure(path.join(outdir,"cosmo"+str(sim)+"_vel_DLA_z"+str(snap)))
    plt.clf()
    (vbin, vels1) = hspec.vel_width_hist("Si", 2, HI_cut = None)
    (vbin, vels2) = hspec.vel_width_hist("Si", 2, HI_cut = 10**20.3)
    mm = np.min((np.size(vels2), np.size(vels1)))
    plt.semilogx(vbin[:mm], vels2[:mm]/vels1[:mm], color="black")
    save_figure(path.join(outdir,"cosmo"+str(sim)+"_vel_DLA_rel_z"+str(snap)))
    plt.clf()

def test_spec_resolution():
    """Plot the velocity widths for different spectral resolutions"""
    #Do spectral resolution test
    halo = myname.get_name(0, False)
    #Higher resolution spectrum
    hspec = ps.PlottingSpectra(5, halo, None, None, savefile="spectra.hdf5")
    hspec.plot_vel_width("Si",2, color="red")
    hspec2 = ps.PlottingSpectra(5, halo, None, None, savefile="spectra4096.hdf5")
    hspec2.plot_vel_width("Si", 2, color="blue")
    vel_data.plot_prochaska_2008_data()
    plt.xlim(1, 1000)
    plt.ylim(1e-5, 2e-2)
    save_figure(path.join(outdir,"cosmo_vel_width_z3_spectra_pix"))
    plt.clf()

    metal_col_den = np.max(hspec.get_col_density("Si", 2),axis=1)
    vel= hspec.vel_width(hspec.metals[("Si",2)][3])
    plt.loglog(metal_col_den, vel, 'o')
    save_figure(path.join(outdir,"cosmo0_vel_col_spectra_pix_low"))
    plt.clf()

    metal_col_den = np.max(hspec2.get_col_density("Si", 2),axis=1)
    vel= hspec.vel_width(hspec2.metals[("Si",2)][3])
    plt.loglog(metal_col_den, vel, 'o')
    save_figure(path.join(outdir,"cosmo0_vel_col_spectra_pix"))
    plt.clf()


#test_spec_resolution()

plot_vel_width_SiII(0, 3)

for ii in (0,3):
    #Plot effect of ignoring low column density metals
    plot_vel_width_metcol(ii, 3)

for ii in (0,3):
    #Plot effect of HI cut
    plot_vel_width_DLA(ii,3)
