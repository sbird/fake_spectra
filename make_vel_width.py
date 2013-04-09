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

def plot_vel_width_sim(sim, snap, color="red", ff=False, HI_cut = None, met_cut = 1e13, unres = 20):
    """Load a simulation and plot its velocity width"""
    halo = "Cosmo"+str(sim)+"_V6"
    if ff:
        halo+="_512"
    #Load from a save file only
    hspec = ps.PlottingSpectra(snap, base+halo, None, None)
    hspec.plot_vel_width("Si", 2, color=color, HI_cut = HI_cut, met_cut = met_cut, unres = unres)

def plot_rel_vel_width(sim1, sim2, snap, color="black"):
    """Load and make a plot of the difference between two simulations"""
    halo1 = "Cosmo"+str(sim1)+"_V6"
    halo2 = "Cosmo"+str(sim2)+"_V6"
    hspec1 = ps.PlottingSpectra(snap, base+halo1, None, None)
    (vbin, vels1) = hspec1.vel_width_hist("Si", 2)
    del hspec1
    hspec1 = ps.PlottingSpectra(snap, base+halo2, None, None)
    (vbin, vels2) = hspec1.vel_width_hist("Si", 2)
    del hspec1
    mm = np.min((np.size(vels2), np.size(vels1)))
    plt.semilogx(vbin[:mm], vels2[:mm]/vels1[:mm], color=color)

def plot_rel_vel_width_DLA(sim, snap, color="red"):
    """Load a simulation and plot its velocity width"""
    halo = "Cosmo"+str(sim)+"_V6"
    #Load from a save file only
    hspec = ps.PlottingSpectra(snap, base+halo, None, None)
    (vbin, vels1) = hspec.vel_width_hist("Si", 2)
    (vbin, vels_d) = hspec.vel_width_hist("Si", 2, HI_cut = 10**20.3)
    mm = np.min((np.size(vels_d), np.size(vels1)))
    plt.semilogx(vbin[:mm], vels_d[:mm]/vels1[:mm], color=color)

def plot_max_col_den(sim, snap, ff=False):
    """Load a simulation and plot the metal column density vs the HI column density"""
    halo = "Cosmo"+str(sim)+"_V6"
    if ff:
        halo+="_512"
    #Load from a save file only
    hspec = ps.PlottingSpectra(snap, base+halo, None, None)
    metal_col_den = np.max(hspec.get_col_density("Si", 2),axis=1)
    HI_col_den = np.max(hspec.get_col_density("H", 1),axis=1)
    plt.loglog(HI_col_den, metal_col_den, 'o')

def plot_spectrum(sim, snap, num):
    """Plot a spectrum"""
    halo = "Cosmo"+str(sim)+"_V6"
    #Load from a save file only
    hspec = ps.PlottingSpectra(snap, base+halo, None, None)
    tau = hspec.get_tau("Si", 2)
    hspec.plot_spectrum(tau, num)

def plot_vel_width_metcol(sim, snap, ff=False):
    """Load a simulation and plot its velocity width"""
    halo = "Cosmo"+str(sim)+"_V6"
    if ff:
        halo+="_512"
    #Load from a save file only
    hspec = ps.PlottingSpectra(snap, base+halo, None, None)
    hspec.plot_vel_width("Si", 2, met_cut = None)
    hspec.plot_vel_width("Si", 2, met_cut = 1e13, color="blue")
    vel_data.plot_prochaska_2008_data()

def plot_vel_width_DLA(sim, snap, ff=False):
    """Plot the effect of the HI cut"""
    halo = "Cosmo"+str(sim)+"_V6"
    #Load from a save file only
    hspec = ps.PlottingSpectra(snap, base+halo, None, None)
    hspec.plot_vel_width("Si", 2, color="red", HI_cut = None)
    hspec.plot_vel_width("Si", 2, color="blue", HI_cut = 10**20.3)
    vel_data.plot_prochaska_2008_data()

def plot_vel_col_den(sim, snap, ff=False):
    """Load a simulation and plot the metal column density vs the HI column density"""
    halo = "Cosmo"+str(sim)+"_V6"
    if ff:
        halo+="_512"
    #Load from a save file only
    hspec = ps.PlottingSpectra(snap, base+halo, None, None)
    metal_col_den = np.max(hspec.get_col_density("Si", 2),axis=1)
    vel= hspec.vel_width(hspec.metals[("Si",2)][3])
    plt.loglog(metal_col_den, vel, 'o')

def test_spec_resolution():
    """Plot the velocity widths for different spectral resolutions"""
    #Do spectral resolution test
    halo = "Cosmo0_V6"
    #Higher resolution spectrum
    hspec = ps.PlottingSpectra(68, base+halo, None, None, savefile=base+halo+"/snapdir_068/spectra.hdf5")
    hspec.plot_vel_width("Si",2, color="red")
    hspec2 = ps.PlottingSpectra(68, base+halo, None, None, savefile=base+halo+"/snapdir_068/spectra4096.hdf5")
    hspec2.plot_vel_width("Si", 2, color="blue", unres=10)
    del hspec2
    vel_data.plot_prochaska_2008_data()
    plt.xlim(10, 1000)
    plt.ylim(1e-5, 2e-2)
    save_figure(path.join(outdir,"cosmo_vel_width_z3_spectra_pix"))
    plt.clf()

    metal_col_den = np.max(hspec.get_col_density("Si", 2),axis=1)
    vel= hspec.vel_width(hspec.metals[("Si",2)][3])
    plt.loglog(metal_col_den, vel, 'o')
    del hspec

colors=["red", "blue", "orange", "black"]

plot_spectrum(0,60, 102)
plt.xlim(1100, 1600)
save_figure(path.join(outdir,"cosmo0_Si_spectrum"))
plt.clf()

#Best-fit base model
plot_vel_width_sim(0, 60, "blue", HI_cut = 10**20.3)
vel_data.plot_prochaska_2008_data()
plt.xlim(10, 1000)
plt.ylim(1e-5, 2e-2)
save_figure(path.join(outdir,"cosmo_vel_width_z3"))
plt.clf()

test_spec_resolution()
save_figure(path.join(outdir,"cosmo0_vel_col_spectra_pix"))
plt.clf()

#The vel widths for different simulations
for ss in (3,0):
    plot_vel_width_sim(ss, 60, colors[ss])

vel_data.plot_prochaska_2008_data()
plt.xlim(10, 1000)
save_figure(path.join(outdir,"cosmo_feedback_z3"))
plt.clf()

#A plot of the redshift evolution
zz = [54,60,68]
for ii in (0,1,2):
    plot_vel_width_sim(0, zz[ii], color=colors[ii])

vel_data.plot_prochaska_2008_data()
plt.xlim(10, 1000)
plt.title("Velocity Widths at z=4-2")
save_figure(path.join(outdir,"cosmo_vel_width_zz"))
plt.clf()

plot_rel_vel_width(1,0,60)
plot_rel_vel_width(2,0,60, color="grey")
#Something odd about this one...
#plot_rel_vel_width(3,0,60, color="red")

save_figure(path.join(outdir,"cosmo_rel_vel_width_z3"))
plt.clf()

for ii in (0,3):
    #Plot effect of ignoring low column density metals
    plot_vel_width_metcol(ii, 60)
    save_figure(path.join(outdir,"cosmo"+str(ii)+"_low_metals_z3"))
    plt.clf()

for ii in (0,3):
    #Plot effect of HI cut
    plot_vel_width_DLA(ii,60)
    save_figure(path.join(outdir,"cosmo"+str(ii)+"_vel_DLA"))
    plt.clf()

for ii in (0,1,2,3):
    #Plot col_density of metals vs HI
    plot_max_col_den(ii, 60)
    save_figure(path.join(outdir,"cosmo"+str(ii)+"z3_coldens"))
    plt.clf()

for ii in (0,1,2,3):
    #Plot metal col. den vs vel width
    plot_vel_col_den(ii, 60)
    save_figure(path.join(outdir,"cosmo"+str(ii)+"_vel_col_z3"))
    plt.clf()

