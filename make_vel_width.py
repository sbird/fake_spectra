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
import sys
from save_figure import save_figure

base="/home/spb/scratch/Cosmo/"
outdir = base + "plots/"
print "Plots at: ",outdir

def plot_vel_width_sim(sim, snap, color="red", ff=False):
    """Load a simulation and plot its velocity width"""
    halo = "Cosmo"+str(sim)+"_V6"
    if ff:
        halo+="_512"
    #Load from a save file only
    hspec = ps.PlottingSpectra(snap, base+halo, None, None)
    tau = hspec.metals[("Si",2)][3]
    metal_col_den = np.max(hspec.get_col_density("Si", 2),axis=1)
    hspec.plot_vel_width(tau, color=color)
    del hspec

def plot_rel_vel_width(sim1, sim2, snap, color="black"):
    """Load and make a plot of the difference between two simulations"""
    halo1 = "Cosmo"+str(sim1)+"_V6"
    halo2 = "Cosmo"+str(sim2)+"_V6"
    hspec1 = ps.PlottingSpectra(snap, base+halo1, None, None)
    (bin, vels1) = hspec1.vel_width_hist(hspec1.metals[("Si",2)][3], 0.1, None)
    del hspec1
    hspec1 = ps.PlottingSpectra(snap, base+halo2, None, None)
    (bin, vels2) = hspec1.vel_width_hist(hspec1.metals[("Si",2)][3], 0.1, None)
    del hspec1
    mm = np.min((np.size(vels2), np.size(vels1)))
    plt.semilogx(bin[:mm], vels2[:mm]/vels1[:mm], color=color)

def plot_rel_vel_width_DLA(sim, snap, color="red"):
    """Load a simulation and plot its velocity width"""
    halo = "Cosmo"+str(sim)+"_V6"
    #Load from a save file only
    hspec = ps.PlottingSpectra(snap, base+halo, None, None)
    (bin, vels1) = hspec.vel_width_hist(hspec.metals[("Si",2)][3], 0.1, None)
    (bin, vels_d) = hspec.vel_width_hist(hspec.metals[("Si",2)][3], 0.1, hspec.get_col_density("H",1))
    mm = np.min((np.size(vels_d), np.size(vels1)))
    plt.semilogx(bin[:mm], vels_d[:mm]/vels1[:mm], color=color)

def plot_vel_width_DLA(sim, snap, color="red", ff=False):
    """Load a simulation and plot its velocity width"""
    halo = "Cosmo"+str(sim)+"_V6"
    if ff:
        halo+="_512"
    #Load from a save file only
    hspec = ps.PlottingSpectra(snap, base+halo, None, None)
    hspec.plot_vel_width(hspec.metals[("Si",2)][3])
    tau = hspec.metals[("Si",2)][3]
    metal_col_den = np.max(hspec.get_col_density("Si", 2),axis=1)
    ind = np.where (metal_col_den > 1e13)
    tau = tau[ind]
    hspec.plot_vel_width(tau, col_rho = hspec.get_col_density("H",1)[ind], color="blue")

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

def plot_vel_width_metcol(sim, snap, ff=False):
    """Load a simulation and plot its velocity width"""
    halo = "Cosmo"+str(sim)+"_V6"
    if ff:
        halo+="_512"
    #Load from a save file only
    hspec = ps.PlottingSpectra(snap, base+halo, None, None)
    hspec.plot_vel_width(hspec.metals[("Si",2)][3])
    tau = hspec.metals[("Si",2)][3]
    metal_col_den = np.max(hspec.get_col_density("Si", 2),axis=1)
    ind = np.where (metal_col_den > 1e13)
    hspec.plot_vel_width(tau[ind],color="blue")

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

colors=["red", "blue", "orange", "purple"]

#The vel widths for different simulations
for ss in (3,2,1,0):
    plot_vel_width_sim(ss, 60, colors[ss])

vel_data.plot_prochaska_2008_data()
plt.xlim(0, 1000)
plt.title("Velocity Widths at z=3")
save_figure(path.join(outdir,"cosmo_vel_width_z3"))
plt.clf()

#Do spectral resolution test
plot_vel_width_sim(0, 60, "red")
halo = "Cosmo0_V6"
#Higher resolution spectrum
hspec = ps.PlottingSpectra(60, base+halo, None, None, savefile=base+halo+"/snapdir_060/spectra2048.hdf5")
tau = hspec.metals[("Si",2)][3]
hspec.plot_vel_width(tau, color="blue")
vel_data.plot_prochaska_2008_data()
plt.xlim(0, 1000)
plt.title("Velocity Widths at z=3")
save_figure(path.join(outdir,"cosmo_vel_width_z3_spectra_pix"))
plt.clf()

metal_col_den = np.max(hspec.get_col_density("Si", 2),axis=1)
vel= hspec.vel_width(hspec.metals[("Si",2)][3])
plt.loglog(metal_col_den, vel, 'o')
save_figure(path.join(outdir,"cosmo0_vel_col_spectra_pix"))
plt.clf()
del hspec

#A plot of the redshift evolution
zz = [54,60,68]
for ii in (0,1,2):
    plot_vel_width_sim(0, zz[ii], color=colors[ii])

vel_data.plot_prochaska_2008_data()
plt.xlim(0, 1000)
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
    plot_vel_width_metcol(ii, 60)
    vel_data.plot_prochaska_2008_data()
    save_figure(path.join(outdir,"cosmo"+str(ii)+"_vel_col_z3"))
    plt.clf()

for i in (0,1,2):
    plot_vel_width_DLA(0,zz[i],ff=True)
    vel_data.plot_prochaska_2008_data()
    save_figure(path.join(outdir,"cosmo512_vel_widthz"+str(zz[i])))
    plt.clf()

for ii in (0,3):
    plot_vel_width_DLA(ii, 60)
    vel_data.plot_prochaska_2008_data()
    save_figure(path.join(outdir,"cosmo"+str(ii)+"_rel_vel_width_DLA"))
    plt.clf()

plot_max_col_den(0, 60,ff=True)
save_figure(path.join(outdir,"cosmo0_512z3_coldens"))
plt.clf()

plot_vel_col_den(0, 60,ff=True)
save_figure(path.join(outdir,"cosmo0_512z3_vvcoldens"))
plt.clf()

plot_vel_width_metcol(0, 60, True)
vel_data.plot_prochaska_2008_data()
save_figure(path.join(outdir,"cosmo0_512_rel_vel_width_metcol"))
plt.clf()

for ii in (0,1,2,3):
    plot_max_col_den(ii, 60)
    save_figure(path.join(outdir,"cosmo"+str(ii)+"z3_coldens"))
    plt.clf()

