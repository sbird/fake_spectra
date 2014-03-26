# -*- coding: utf-8 -*-
"""Make some plots of 2d histograms, velocity widths vs some other things"""

import matplotlib
matplotlib.use('PDF')

import matplotlib.pyplot as plt

import plot_spectra as ps
import os.path as path
import myname
import numpy as np
from save_figure import save_figure

outdir = path.join(myname.base, "plots/mass/")
print "Plots at: ",outdir

zrange = {1:(7,3.5), 3:(7,0), 5:(2.5,0)}
#Colors and linestyles for the simulations
colors = {0:"red", 1:"purple", 2:"cyan", 3:"green", 4:"gold", 5:"orange", 7:"blue", 6:"grey"}
colors2 = {0:"darkred", 1:"indigo", 2:"cyan", 3:"darkgreen", 4:"gold", 5:"orange", 7:"darkblue", 6:"grey"}
lss = {0:"--",1:":", 2:":",3:"-.", 4:"--", 5:"-",6:"--",7:"-"}
labels = {0:"REF",1:"HVEL", 2:"HVNA",3:"NOSN", 4:"NAWW", 5:"MVEL",6:"METAL",7:"2xUV"}

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
    plt.xlim(10,300)

def plot_vir_vsmass(sim, snap):
    """Plot mass histogram"""
    #Load from a save file only
    hspec = get_hspec(sim,snap)
    (halos, _) = hspec.find_nearest_halo()
    hspec._plot_xx_vs_mass(hspec.sub_mass[halos], name = "Mass", color=colors[sim], color2=colors2[sim], log=True)
    out = "cosmo"+str(sim)+"_vir_mass_z"+str(snap)
    save_figure(path.join(outdir,out))
    plt.clf()

def plot_vvir(sim, snap):
    """Plot histogram of velocity width by virial velocity"""
    #Load from a save file only
    hspec = get_hspec(sim,snap)
    (mbins, pdf) = hspec.mass_hist()
    hspec.plot_virial_vel_vs_vel_width("Si", 2, color=colors[sim], ls=lss[sim], label=labels[sim])

def plot_mass_vs(sim, snap):
    """Plot mass vs metallicity and vel width"""
    out = "cosmo"+str(sim)+"_met_mass_z"+str(snap)
    #Load from a save file only
    hspec = get_hspec(sim,snap)
    hspec.plot_Z_vs_mass(color=colors[sim], color2=colors2[sim])
    plt.ylabel(r"$Z (Z_\odot) $")
    plt.xlabel(r"$v_\mathrm{vir}$")
    save_figure(path.join(outdir,out))
    plt.clf()
    out = "cosmo"+str(sim)+"_vel_mass_z"+str(snap)
    hspec.plot_vel_vs_mass("Si",2, color=colors[sim], color2=colors2[sim])
    plt.loglog(np.logspace(0.5,3.5), np.logspace(0.5,3.5), ls="-", color="black")
    plt.ylabel(r"$v_\mathrm{90}$ (km s$^{-1}$)")
    plt.xlabel(r"$v_\mathrm{vir}$ (km s$^{-1}$)")
    plt.xlim(10,300)
    plt.ylim(10,500)
    save_figure(path.join(outdir,out))
    plt.clf()

def plot_mass_vs_mm(sim, snap):
    """Plot mass vs mm and fedge"""
    out = "cosmo"+str(sim)+"_fmm_mass_z"+str(snap)
    #Load from a save file only
    hspec = get_hspec(sim,snap)
    fmm = hspec.vel_mean_median("Si",2)
    hspec._plot_xx_vs_mass(fmm, name = "fmm", color=colors[sim], color2=colors2[sim], log=False)
    plt.ylabel(r"$f_\mathrm{mm}$")
    plt.xlabel(r"$v_\mathrm{vir}$ (km s$^{-1}$)")
    save_figure(path.join(outdir,out))
    plt.clf()
    out = "cosmo"+str(sim)+"_fedge_mass_z"+str(snap)
    fpk = hspec.vel_peak("Si",2)
    hspec._plot_xx_vs_mass(fpk, name = "fedge", color=colors[sim], color2=colors2[sim], log=False)
    plt.ylabel(r"$f_\mathrm{edge}$")
    plt.xlabel(r"$v_\mathrm{vir}$ (km s$^{-1}$)")
    plt.xlim(10,300)
    save_figure(path.join(outdir,out))
    plt.clf()
    out = "cosmo"+str(sim)+"_eqw_mass_z"+str(snap)
    eq_width = hspec.equivalent_width("Si", 2, 1526)
    hspec._plot_xx_vs_mass(eq_width, name = "eqw", color=colors[sim], color2=colors2[sim])
    plt.ylabel(r"$W_\mathrm{1526}$")
    plt.xlabel(r"$v_\mathrm{vir}$ (km s$^{-1}$)")
    plt.xlim(10,300)
    save_figure(path.join(outdir,out))
    plt.clf()

def plot_mm_vs_vel(sim, snap):
    """Plot vel width vs mm and fedge"""
    hspec = get_hspec(sim,snap)
    out = "cosmo"+str(sim)+"_fmm_vel_z"+str(snap)
    vel = hspec.vel_width("Si", 2)
    ii = hspec.get_filt("Si", 2)
    fmm = hspec.vel_mean_median("Si",2)
    fpk = hspec.vel_peak("Si",2)
#     (halo, _) = hspec.find_nearest_halo()
#     ind2 = np.where(halo > 0)
#     halo = halo[ind2]
#     vel = vel[ind2]
#     fpk = fpk[ind2]
#     fmm = fmm[ind2]
#     virial = hspec.virial_vel(halo)
    hspec._plot_2d_contour(vel[ii], fmm[ii], 10, "fmm vel", color=colors[sim], color2=colors2[sim], ylog=False)
    plt.xlim(10,2e3)
    save_figure(path.join(outdir,out))
    plt.clf()
    out = "cosmo"+str(sim)+"_fedge_vel_z"+str(snap)
    hspec._plot_2d_contour(vel[ii], fpk[ii], 10, "fpk vel", color=colors[sim], color2=colors2[sim], ylog=False)
    plt.ylabel(r"$f_\mathrm{edge}$")
    plt.xlabel(r"$v_\mathrm{90}$ (km s$^{-1}$)")
    plt.xlim(10,300)
    save_figure(path.join(outdir,out))
    plt.clf()
#     out = "cosmo"+str(sim)+"_fedge_fmm_z"+str(snap)
#     hspec._plot_2d_contour(fmm, fpk, 10, "fpk fmm", color=colors[sim], color2=colors2[sim], ylog=False, xlog=False)
#     save_figure(path.join(outdir,out))
#     plt.clf()

def plot_mult_frac(sim, snap):
    """Plot fraction of lines from separated halos"""
    hspec = get_hspec(sim, snap)
    hspec.plot_mult_halo_frac(color=colors[sim], color2 = colors2[sim], ls=lss[sim])

if __name__ == "__main__":

    simlist = (0,1,3,7)  #range(8)
    for zz in (1,3,5):
        for ss in simlist:
            plot_mass_hists(ss, zz)
        save_figure(path.join(outdir,"cosmo_halos_feedback_z"+str(zz)))
        plt.clf()

    for zz in (1,3,5):
        for ss in simlist:
            plot_mult_frac(ss,zz)
        plt.legend(loc=2,ncol=3)
        save_figure(path.join(outdir,"cosmo_mult_frac_z"+str(zz)))
        plt.clf()

    for zz in (1,3,5):
        for ss in simlist:
            plot_vvir(ss, zz)
        plt.xlabel(r"$v_\mathrm{90} / v_\mathrm{vir}$")
        plt.xlim(0.05, 30)
        save_figure(path.join(outdir,"cosmo_vw_vel_vir_z"+str(zz)))
        plt.clf()

    for zz in (1,3,5):
        for ss in simlist:
#             plot_vir_vsmass(ss,zz)
            plot_mass_vs(ss, zz)
            plot_mass_vs_mm(ss, zz)
            plot_mm_vs_vel(ss, zz)

