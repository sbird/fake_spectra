# -*- coding: utf-8 -*-
"""Make some plots of 2d histograms, velocity widths vs some other things"""

import matplotlib
matplotlib.use('PDF')

import matplotlib.pyplot as plt

import plot_spectra as ps
import os.path as path
import myname
from save_figure import save_figure

outdir = path.join(myname.base, "plots/mass/")
print "Plots at: ",outdir

zrange = {1:(7,3.5), 3:(7,0), 5:(2.5,0)}
#Colors and linestyles for the simulations
colors = {0:"red", 1:"purple", 2:"cyan", 3:"green", 4:"gold", 5:"orange", 7:"blue", 6:"grey"}
lss = {0:"--",1:":", 2:":",3:"-.", 4:"--", 5:"-",6:"--",7:"-"}
labels = {0:"REF",1:"HVEL", 2:"HVNA",3:"NOSN", 4:"NAWW", 5:"MVEL",6:"METAL",7:"TUV"}

hspec_cache = {}

def get_hspec(sim, snap):
    """Get a spectra object, possibly from the cache"""
    halo = myname.get_name(sim, True)
    #Load from a save file only
    try:
        hspec = hspec_cache[(sim, snap)]
    except KeyError:
        hspec = ps.PlottingSpectra(snap, halo, label=labels[sim])
        hspec_cache[(sim, snap)] = hspec
    return hspec


def plot_mass_hists(sim, snap):
    """Plot mass histogram"""
    #Load from a save file only
    hspec = get_hspec(sim,snap)
    (mbins, pdf) = hspec.mass_hist()
    plt.semilogx(mbins,pdf,color=colors[sim], ls=lss[sim],label=labels[sim])

def plot_mass_vs(sim, snap):
    """Plot mass vs metallicity and vel width"""
    out = "cosmo"+str(sim)+"_met_mass_z"+str(snap)
    #Load from a save file only
    hspec = get_hspec(sim,snap)
    hspec.plot_Z_vs_mass(color=colors[sim])
    save_figure(path.join(outdir,out))
    plt.clf()
    out = "cosmo"+str(sim)+"_vel_mass_z"+str(snap)
    hspec.plot_vel_vs_mass("Si",2, color=colors[sim])
    save_figure(path.join(outdir,out))
    plt.clf()

if __name__ == "__main__":
    for ss in range(8):
        plot_mass_hists(ss, 3)
    save_figure(path.join(outdir,"cosmo_halos_feedback_z3"))
    plt.clf()

    for zz in (1,3,5):
        for ss in range(8):
            plot_mass_vs(ss, zz)

