#!/usr/bin env python
# -*- coding: utf-8 -*-
"""Make some plots of 2d histograms, velocity widths vs some other things"""

import matplotlib
matplotlib.use('PDF')

import matplotlib.pyplot as plt

import vw_plotspectra as ps
import os.path as path
import numpy as np
from save_figure import save_figure
import myname

outdir = path.join(myname.base, "plots/2d_hist/")
print "Plots at: ",outdir

def plot_max_den(sim, snap, ff=True):
    """Load a simulation and plot the max metal density vs the max HI density"""
    halo = myname.get_name(sim,ff)
    #Load from a save file only
    hspec = ps.VWPlotSpectra(snap, halo, None, None)
    metal_den = np.max(hspec.get_density("Si", 2),axis=1)
    HI_den = np.max(hspec.get_density("H", 1),axis=1)
    ind = hspec.get_filt("Si",2)
    (H, xedges, yedges) = np.histogram2d(np.log10(metal_den[ind]), np.log10(HI_den[ind]), bins=30,normed=True)
    extent = [yedges[0], yedges[-1], xedges[-1], xedges[0]]
    plt.imshow(H, extent=extent, aspect="auto", vmax = 0.15)
    plt.colorbar()

def plot_vel_den(sim, snap, ff=True):
    """Load a simulation and plot the metal density vs the velocity width"""
    halo = myname.get_name(sim, ff)
    #Load from a save file only
    hspec = ps.VWPlotSpectra(snap, halo, None, None)
    vel = hspec.vel_width("Si",2)
    den = hspec.get_density("Si", 2)
    ind = hspec.get_filt("Si",2)
    (H, xedges, yedges) = np.histogram2d(np.log10(den[ind]), np.log10(vel[ind]), bins=30,normed=True)
    extent = [yedges[0], yedges[-1], xedges[-1], xedges[0]]
    plt.imshow(H, extent=extent, aspect="auto")
    plt.colorbar()

def plot_vel_HI_col_den(sim, snap, ff=True):
    """Load a simulation and plot the HI column density vs the velocity width"""
    halo = myname.get_name(sim, ff)
    #Load from a save file only
    hspec = ps.VWPlotSpectra(snap, halo, None, None)
    HI_den = np.sum(hspec.get_col_density("H", 1),axis=1)
    vel= hspec.vel_width("Si",2)
    ind = hspec.get_filt("Si",2)
    (H, xedges, yedges) = np.histogram2d(np.log10(HI_den[ind]), np.log10(vel[ind]), bins=30,normed=True)
    extent = [yedges[0], yedges[-1], xedges[-1], xedges[0]]
    plt.imshow(H, extent=extent, aspect="auto")
    plt.colorbar()

def plot_vel_mass(sim, snap, ff=True):
    """Load a simulation and plot the halo mass vs the velocity width"""
    halo = myname.get_name(sim, ff)
    #Load from a save file only
    hspec = ps.VWPlotSpectra(snap, halo)
    ind = hspec.get_filt("Si",2)
    vel= hspec.vel_width("Si",2)[ind]
    (halos, _) = hspec.find_nearest_halo()
    mass = hspec.sub_mass[halos][ind]
    (H, xedges, yedges) = np.histogram2d(np.log10(mass), np.log10(vel), bins=30,normed=True)
    extent = [yedges[0], yedges[-1], xedges[-1], xedges[0]]
    plt.imshow(H, extent=extent, aspect="auto")
    plt.colorbar()

def plot_met_mass(sim, snap, ff=True):
    """Load a simulation and plot the halo mass vs the velocity width"""
    halo = myname.get_name(sim, ff)
    #Load from a save file only
    hspec = ps.VWPlotSpectra(snap, halo)
    ind = hspec.get_filt("Si",2)
    met = hspec.get_metallicity()[ind]
    (halos, _) = hspec.find_nearest_halo()
    mass = hspec.sub_mass[halos][ind]
    (H, xedges, yedges) = np.histogram2d(np.log10(mass), np.log10(met), bins=30,normed=True)
    extent = [yedges[0], yedges[-1], xedges[-1], xedges[0]]
    plt.imshow(H, extent=extent, aspect="auto")
    plt.colorbar()

def plot_vel_metals(sim, snap, ff=True):
    """Plot the correlation between metallicity and velocity width"""
    halo = myname.get_name(sim, ff)
    #Load from a save file only
    hspec = ps.VWPlotSpectra(snap, halo)
    met = hspec.get_metallicity()
    vel = hspec.vel_width("Si",2)
    (H, xedges, yedges) = np.histogram2d(np.log10(met), np.log10(vel), bins=30,normed=True)
    extent = [yedges[0], yedges[-1], xedges[-1], xedges[0]]
    plt.imshow(H, extent=extent, aspect="auto")
    plt.colorbar()

def plot_Si_metals(sim, snap, ff=True):
    """Plot the correlation between metallicity and velocity width"""
    halo = myname.get_name(sim, ff)
    #Load from a save file only
    hspec = ps.VWPlotSpectra(snap, halo)
    met = hspec.get_metallicity()
    MM = hspec.get_density("Si",2)
    HH = hspec.get_density("H",-1)
    mms = np.sum(MM, axis=1)
    hhs = np.sum(HH, axis=1)
    Simet = mms/hhs/0.0133
    (H, xedges, yedges) = np.histogram2d(np.log10(met), np.log10(Simet), bins=30,normed=True)
    extent = [yedges[0], yedges[-1], xedges[-1], xedges[0]]
    plt.imshow(H, extent=extent, aspect="auto")
    plt.colorbar()

reds = {1:4, 3:3, 5:2}

for ii in (0,1,2,3):
    #Plot col_density of metals vs HI
    plot_vel_mass(ii, 3)
    save_figure(path.join(outdir,"cosmo"+str(ii)+"z3_vel_mass"))
    plt.clf()

for ii in (0,1,2,3):
    #Plot col_density of metals vs HI
    plot_met_mass(ii, 3)
    save_figure(path.join(outdir,"cosmo"+str(ii)+"z3_met_mass"))
    plt.clf()

# plot_Si_metals(0, 3)
# save_figure(path.join(outdir,"cosmo0_512_z3_Si_metals"))
# plt.clf()
#

for ii in (0,1,2,3):
    #Plot col_density of metals vs HI
    plot_max_den(ii, 3)
    save_figure(path.join(outdir,"cosmo"+str(ii)+"z3_coldens"))
    plt.clf()

for ii in (0,1,2,3):
    #Plot metal col. den vs vel width
    plot_vel_den(ii, 3)
    save_figure(path.join(outdir,"cosmo"+str(ii)+"_vel_den_z3"))
    plt.clf()

for ii in (0,1,2,3):
    #Plot metal col. den vs vel width
    plot_vel_HI_col_den(ii, 3)
    save_figure(path.join(outdir,"cosmo"+str(ii)+"_vel_HI_col_z3"))
    plt.clf()
