# -*- coding: utf-8 -*-
"""Make some HI related plots from the cosmo runs"""

import matplotlib
matplotlib.use('PDF')

import matplotlib.pyplot as plt

import plot_spectra as ps
import dla_data
import os.path as path
import myname
from save_figure import save_figure

outdir = path.join(myname.base,"plots/spectra_HI")
print "Plots at: ",outdir

def plot_cddf_a_halo(sim, snap, color="red", ff=True):
    """Load a simulation and plot its cddf"""
    halo = myname.get_name(sim, ff)
    hspec = ps.PlottingSpectra(snap, halo)
    hspec.plot_cddf(color=color)
    del hspec


def plot_Omega_DLA(sim, color="red", ff=True):
    """Plot Omega_DLA over a range of redshifts"""
    halo = myname.get_name(sim, ff)
    om = {}
    for snap in (1,3,5):
        hspec = ps.PlottingSpectra(snap, halo)
        om[hspec.red] = hspec.omega_DLA()
    plt.semilogy(om.keys(), om.values(), 'o-', color=color)
    plt.xlabel("z")
    plt.ylabel(r"$\Omega_{DLA}$")
    return om

def plot_rho_HI(sim, color="red", ff=True):
    """Plot rho_HI across redshift"""
    halo = myname.get_name(sim, ff)
    zzz = {4:1, 3:3, 2:5}
    rho_HI = {}
    for zz in (4,3,2):
        try:
            hspec = ps.PlottingSpectra(zzz[zz], halo)
            rho_HI[zz]=hspec.omega_DLA()
            del hspec
        except TypeError:
            pass
    plt.plot(rho_HI.keys(),rho_HI.values(), color=color)

def plot_dndx(sim, color="red", ff=True):
    """Plot dndx (cross-section) across redshift"""
    halo = myname.get_name(sim, ff)
    zzz = {4:1, 3:3, 2:5}
    dndx={}
    for zz in (4,3,2):
        try:
            hspec = ps.PlottingSpectra(zzz[zz], halo)
            dndx[zz]=hspec.line_density()
            del hspec
        except TypeError:
            pass
    plt.plot(dndx.keys(),dndx.values(), color=color)

colors = {0:"red", 1:"purple", 2:"blue", 3:"green", 4:"orange"}

for i in (0,1,2,3,4):
    plot_dndx(i,colors[i])
dla_data.dndx_not()
save_figure(path.join(outdir,"cosmo_dndx"))
plt.clf()

for i in (0,1,2,3,4):
    plot_rho_HI(i,colors[i])
dla_data.omegahi_not()
save_figure(path.join(outdir,"cosmo_rhohi"))
plt.clf()

#Make a plot of the column density functions.
for ss in (4,3,2,1,0):
    plot_cddf_a_halo(ss, 3, color=colors[ss])

dla_data.column_density_data(moment=True)

save_figure(path.join(outdir,"cosmo_cddf_z3"))
plt.clf()

