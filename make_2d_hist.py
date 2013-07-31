#!/usr/bin env python
# -*- coding: utf-8 -*-
"""Make some plots of 2d histograms, velocity widths vs some other things"""

import matplotlib
matplotlib.use('PDF')

import matplotlib.pyplot as plt

import plot_spectra as ps
import os.path as path
import numpy as np
from save_figure import save_figure

base="/home/spb/scratch/Cosmo/"
outdir = base + "plots/2d_hist/"
print "Plots at: ",outdir

def plot_max_col_den(sim, snap, ff=False):
    """Load a simulation and plot the metal column density vs the HI column density"""
    halo = "Cosmo"+str(sim)+"_V6"
    if ff:
        halo+="_512"
    #Load from a save file only
    hspec = ps.PlottingSpectra(snap, base+halo, None, None)
    metal_col_den = np.max(hspec.get_col_density("Si", 2),axis=1)
    HI_col_den = np.max(hspec.get_col_density("H", 1),axis=1)
    ind = np.where(metal_col_den > 1e12)
    (H, xedges, yedges) = np.histogram2d(np.log10(metal_col_den[ind]), np.log10(HI_col_den[ind]), bins=30,normed=True)
    extent = [yedges[0], yedges[-1], xedges[-1], xedges[0]]
    plt.imshow(H, extent=extent, aspect="auto", vmax = 0.15)
    plt.colorbar()

def plot_vel_col_den(sim, snap, ff=False):
    """Load a simulation and plot the metal column density vs the HI column density"""
    halo = "Cosmo"+str(sim)+"_V6"
    if ff:
        halo+="_512"
    #Load from a save file only
    hspec = ps.PlottingSpectra(snap, base+halo, None, None)
    metal_col_den = np.max(hspec.get_col_density("Si", 2),axis=1)
    vel= hspec.vel_width(hspec.get_observer_tau("Si",2))
    ind = np.where(metal_col_den > 1e12)
    (H, xedges, yedges) = np.histogram2d(np.log10(metal_col_den[ind]), np.log10(vel[ind]), bins=30,normed=True)
    extent = [yedges[0], yedges[-1], xedges[-1], xedges[0]]
    plt.imshow(H, extent=extent, aspect="auto")
    plt.colorbar()

def plot_vel_den(sim, snap, ff=False):
    """Load a simulation and plot the metal column density vs the HI column density"""
    halo = "Cosmo"+str(sim)+"_V6"
    if ff:
        halo+="_512"
    #Load from a save file only
    hspec = ps.PlottingSpectra(snap, base+halo, None, None)
    vel = hspec.vel_width(hspec.get_observer_tau("Si",2))
    den = hspec.vel_width(hspec.get_col_density("Si", 2))
    ind = hspec.get_filt("Si",2)
    (H, xedges, yedges) = np.histogram2d(np.log10(den[ind]), np.log10(vel[ind]), bins=30,normed=True)
    extent = [yedges[0], yedges[-1], xedges[-1], xedges[0]]
    plt.imshow(H, extent=extent, aspect="auto")
    plt.colorbar()

def plot_vel_HI_col_den(sim, snap, ff=False):
    """Load a simulation and plot the metal column density vs the HI column density"""
    halo = "Cosmo"+str(sim)+"_V6"
    if ff:
        halo+="_512"
    #Load from a save file only
    hspec = ps.PlottingSpectra(snap, base+halo, None, None)
    metal_col_den = np.max(hspec.get_col_density("Si", 2),axis=1)
    HI_col_den = np.max(hspec.get_col_density("H", 1),axis=1)
    vel= hspec.vel_width(hspec.get_observer_tau("Si",2))
    ind = np.where(metal_col_den > 1e12)
    (H, xedges, yedges) = np.histogram2d(np.log10(HI_col_den[ind]), np.log10(vel[ind]), bins=30,normed=True)
    extent = [yedges[0], yedges[-1], xedges[-1], xedges[0]]
    plt.imshow(H, extent=extent, aspect="auto")
    plt.colorbar()

def plot_vel_mass(sim, snap, ff=False):
    """Load a simulation and plot the halo mass vs the velocity width"""
    halo = "Cosmo"+str(sim)+"_V6"
    if ff:
        halo+="_512"
    #Load from a save file only
    hspec = ps.PlottingSpectra(snap, base+halo)
    ind = hspec.get_filt("Si",2)
    vel= hspec.vel_width(hspec.get_observer_tau("Si",2))[ind]
    (halos, dists) = hspec.find_nearest_halo()
    mass = hspec.sub_mass[halos][ind]
    ind = np.where(vel > 10)
    (H, xedges, yedges) = np.histogram2d(np.log10(mass[ind]), np.log10(vel[ind]), bins=30,normed=True)
    extent = [yedges[0], yedges[-1], xedges[-1], xedges[0]]
    plt.imshow(H, extent=extent, aspect="auto")
    plt.colorbar()

def plot_met_mass(sim, snap, ff=False):
    """Load a simulation and plot the halo mass vs the velocity width"""
    halo = "Cosmo"+str(sim)+"_V6"
    if ff:
        halo+="_512"
    #Load from a save file only
    hspec = ps.PlottingSpectra(snap, base+halo)
    ind = hspec.get_filt("Si",2)
    met = hspec.get_metallicity()[ind]
    (halos, dists) = hspec.find_nearest_halo()
    mass = hspec.sub_mass[halos][ind]
    (H, xedges, yedges) = np.histogram2d(np.log10(mass), np.log10(met), bins=30,normed=True)
    extent = [yedges[0], yedges[-1], xedges[-1], xedges[0]]
    plt.imshow(H, extent=extent, aspect="auto")
    plt.colorbar()

def plot_vel_metals(sim, snap, ff=False):
    """Plot the correlation between metallicity and velocity width"""
    halo = "Cosmo"+str(sim)+"_V6"
    if ff:
        halo+="_512"
    #Load from a save file only
    hspec = ps.PlottingSpectra(snap, base+halo)
    met = hspec.get_metallicity()
    tau = hspec.get_observer_tau("Si", 2)
    vel = hspec.vel_width(tau)
    #Ignore objects too faint to be seen or unresolved
    ind2 = np.where(np.logical_and(vel > 15, met > 1e-3))
    (H, xedges, yedges) = np.histogram2d(np.log10(met[ind2]), np.log10(vel[ind2]), bins=30,normed=True)
    extent = [yedges[0], yedges[-1], xedges[-1], xedges[0]]
    plt.imshow(H, extent=extent, aspect="auto")
    plt.colorbar()

def plot_Si_metals(sim, snap, ff=False):
    """Plot the correlation between metallicity and velocity width"""
    halo = "Cosmo"+str(sim)+"_V6"
    if ff:
        halo+="_512"
    #Load from a save file only
    hspec = ps.PlottingSpectra(snap, base+halo)
    met = hspec.get_metallicity()
    MM = hspec.get_col_density("Si",2)
    HH = hspec.get_col_density("H",-1)
    mms = np.sum(MM, axis=1)
    hhs = np.sum(HH, axis=1)
    Simet = mms/hhs/0.0133
    #Ignore objects too faint to be seen or unresolved
    ind2 = np.where(met > 1e-3)
    (H, xedges, yedges) = np.histogram2d(np.log10(met[ind2]), np.log10(Simet[ind2]), bins=30,normed=True)
    extent = [yedges[0], yedges[-1], xedges[-1], xedges[0]]
    plt.imshow(H, extent=extent, aspect="auto")
    plt.colorbar()

reds = {1:4, 3:3, 5:2}

plot_vel_mass(0, 3,True)
save_figure(path.join(outdir,"cosmo0_512_z3_vel_mass"))
plt.clf()

plot_met_mass(0, 3,True)
save_figure(path.join(outdir,"cosmo0_512_z3_met_mass"))
plt.clf()

import sys
sys.exit()
plot_vel_metals(0, 3,True)
save_figure(path.join(outdir,"cosmo0_512_z3_metals"))
plt.clf()

plot_Si_metals(0, 3,True)
save_figure(path.join(outdir,"cosmo0_512_z3_Si_metals"))
plt.clf()


for ii in (0,1,2,3):
    #Plot col_density of metals vs HI
    plot_max_col_den(ii, 3)
    save_figure(path.join(outdir,"cosmo"+str(ii)+"z3_coldens"))
    plt.clf()

for ii in (0,1,2,3):
    #Plot metal col. den vs vel width
    plot_vel_den(ii, 3)
    save_figure(path.join(outdir,"cosmo"+str(ii)+"_vel_den_z3"))
    plt.clf()


for ii in (0,1,2,3):
    #Plot metal col. den vs vel width
    plot_vel_col_den(ii, 3)
    save_figure(path.join(outdir,"cosmo"+str(ii)+"_vel_col_z3"))
    plt.clf()

for ii in (0,1,2,3):
    #Plot metal col. den vs vel width
    plot_vel_HI_col_den(ii, 3)
    save_figure(path.join(outdir,"cosmo"+str(ii)+"_vel_HI_col_z3"))
    plt.clf()
