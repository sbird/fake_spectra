#!/usr/bin env python
# -*- coding: utf-8 -*-
"""Make some HI related plots from the cosmo runs"""

import matplotlib
matplotlib.use('PDF')

import matplotlib.pyplot as plt

import plot_spectra as ps
import dla_data
import os.path as path
from save_figure import save_figure

base="/home/spb/scratch/Cosmo/"
outdir = base + "plots/spectra_HI"
print "Plots at: ",outdir

def plot_cddf_a_halo(sim, snap, color="red", ff=False):
    """Load a simulation and plot its cddf"""
    halo = "Cosmo"+str(sim)+"_V6"
    if ff:
        halo+="_512"
    hspec = ps.PlottingSpectra(snap, base+halo, savefile="rand_spectra.hdf5")
    hspec.plot_cddf(color=color)
    del hspec


def plot_Omega_DLA(sim, color="red", ff=False):
    """Plot Omega_DLA over a range of redshifts"""
    halo = "Cosmo"+str(sim)+"_V6"
    if ff:
        halo+="_512"
    om = {}
    for snap in (54, 60, 68):
        hspec = ps.PlottingSpectra(snap, base+halo, savefile="rand_spectra.hdf5")
        om[hspec.red] = hspec.omega_DLA()
    plt.semilogy(om.keys(), om.values(), 'o-', color=color)
    plt.xlabel("z")
    plt.ylabel(r"$\Omega_{DLA}$")
    return om

def plot_rho_HI(sim, color="red", ff=False):
    """Plot rho_HI across redshift"""
    halo = "Cosmo"+str(sim)+"_V6"
    if ff:
        halo+="_512"
    zzz = {4:54, 3:60, 2:68}
    rho_HI = {}
    for zz in (4,3,2):
        try:
            hspec = ps.PlottingSpectra(zzz[zz], base+halo, savefile="rand_spectra.hdf5")
            rho_HI[zz]=hspec.rho_DLA()
            del hspec
        except TypeError:
            pass
    plt.plot(rho_HI.keys(),rho_HI.values(), color=color)

def plot_dndx(sim, color="red", ff=False):
    """Plot dndx (cross-section) across redshift"""
    halo = "Cosmo"+str(sim)+"_V6"
    if ff:
        halo+="_512"
    zzz = {4:54, 3:60, 2:68}
    dndx={}
    for zz in (4,3,2):
        try:
            hspec = ps.PlottingSpectra(zzz[zz], base+halo, savefile="rand_spectra.hdf5")
            dndx[zz]=hspec.line_density()
            del hspec
        except TypeError:
            pass
    plt.plot(dndx.keys(),dndx.values(), color=color)

colors=["red", "blue", "orange", "purple"]

for i in (0,2,3):
    plot_dndx(i,colors[i])
plot_dndx(0,colors[1], True)
dla_data.dndx()
save_figure(path.join(outdir,"cosmo_dndx"))
plt.clf()

for i in (0,2,3):
    plot_rho_HI(i,colors[i])
plot_rho_HI(0,colors[1], True)
dla_data.rhohi()
save_figure(path.join(outdir,"cosmo_rhohi"))
plt.clf()

#Make a plot of the column density functions.
for ss in (3,2,0):
    plot_cddf_a_halo(ss, 60, color=colors[ss])
plot_cddf_a_halo(0, 60, colors[1],True)

dla_data.column_density_data()

save_figure(path.join(outdir,"cosmo_cddf_z3"))
plt.clf()

