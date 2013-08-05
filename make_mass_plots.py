# -*- coding: utf-8 -*-
"""Make some plots of 2d histograms, velocity widths vs some other things"""

import matplotlib
matplotlib.use('PDF')

import matplotlib.pyplot as plt

import plot_spectra as ps
import os.path as path
import myname
from save_figure import save_figure

colors = {0:"red", 1:"purple", 2:"blue", 3:"green", 4:"orange"}
lss = {0:"--",1:":", 2:"-",3:"-.", 4:"-"}
outdir = path.join(myname.base, "plots/mass/")
print "Plots at: ",outdir


def plot_mass_hists(sim, snap, ff=True):
    """Plot mass histogram"""
    halo = myname.get_name(sim, ff)
    #Load from a save file only
    hspec = ps.PlottingSpectra(snap, halo)
    (mbins, pdf) = hspec.mass_hist()
    plt.semilogx(mbins,pdf,color=colors[sim], ls=lss[sim])

def plot_mass_vs(sim, snap, ff=True):
    """Plot mass vs metallicity and vel width"""
    halo = myname.get_name(sim, ff)
    out = "cosmo"+str(sim)+"_met_mass_z"+str(snap)
    if ff:
        out+="_512"
    #Load from a save file only
    hspec = ps.PlottingSpectra(snap, halo)
    hspec.plot_Z_vs_mass(color=colors[sim])
    save_figure(path.join(outdir,out))
    plt.clf()
    out = "cosmo"+str(sim)+"_vel_mass_z"+str(snap)
    if ff:
        out+="_512"
    hspec.plot_vel_vs_mass("Si",2, color=colors[sim])
    save_figure(path.join(outdir,out))
    plt.clf()

if __name__ == "__main__":
    for ss in (0,1,2,3,4):
        plot_mass_hists(ss, 3)
    save_figure(path.join(outdir,"cosmo_halos_feedback_z3"))
    plt.clf()

    for zz in (1,3,5):
        plot_mass_hists(1,zz)
    save_figure(path.join(outdir,"cosmo_halos_1_zz"))
    plt.clf()

    for ss in (0,1,2,3,4):
        for zz in (1,3,5):
            plot_mass_vs(ss, zz, True)

