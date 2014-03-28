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

zrange = {1:(7,3.5), 3:(7,0), 5:(2.5,0)}

def plot_metal_ion_corr(sim, snap,species="Si",ion=2):
    """Plot metallicity from Z/H vs from a single species for computing ionisation corrections"""
    halo = myname.get_name(sim)
    hspec = ps.PlottingSpectra(snap, halo)
    hspec.plot_metallicity(color="red", ls="-")
    hspec.plot_species_metallicity(species, ion, color="blue", ls="-")
    vel_data.plot_alpha_metal_data((3.5,2.5))
    save_figure(path.join(outdir, "cosmo"+str(sim)+"_ion_corr"+str(snap)))
    plt.clf()
    hspec.plot_ion_corr(species, ion)
    save_figure(path.join(outdir, "cosmo"+str(sim)+"_rel_ion_corr"+str(snap)))
    plt.clf()

def plot_vel_width_metcol(sim, snap):
    """Load a simulation and plot its velocity width"""
    halo = myname.get_name(sim)
    #Load from a save file only
    hspec = ps.PlottingSpectra(snap, halo)
    hspec.plot_vel_width("Si", 2)
    hspec.plot_vel_width("Si", 2, color="blue")
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

def test_spec_resolution():
    """Plot the velocity widths for different spectral resolutions"""
    halo = myname.get_name(7)
    #Higher resolution spectrum
    hspec = ps.PlottingSpectra(3, halo, savefile="grid_spectra_DLA.hdf5")
    hspec2 = ps.PlottingSpectra(3, halo, savefile="grid_spectra_DLA_res.hdf5")
    plot_check(hspec,hspec2,"specres")

def test_pecvel():
    """Plot the velocity widths with and without peculiar velocities"""
    halo = myname.get_name(7)
    #Higher resolution spectrum
    hspec = ps.PlottingSpectra(3, halo, savefile="grid_spectra_DLA.hdf5")
    hspec2 = ps.PlottingSpectra(3, halo, None, None, savefile="grid_spectra_DLA_pecvel.hdf5")
    plot_check(hspec,hspec2,"pecvel")

def test_tophat():
    """Plot the velocity widths with and with top hat vs SPH"""
    halo = myname.get_name(7)
    #Higher resolution spectrum
    hspec = ps.PlottingSpectra(3, halo, savefile="grid_spectra_DLA.hdf5")
    hspec2 = ps.PlottingSpectra(3, halo, None, None, savefile="grid_spectra_DLA_tophat.hdf5")
    plot_check(hspec,hspec2,"tophat")

def test_box_resolution():
    """Plot the velocity widths for different size boxes"""
    halo = myname.get_name(5)
    halo10 = myname.get_name(5,box=10)
    for zz in (1,3,5):
        hspec = ps.PlottingSpectra(zz, halo, label="MVEL")
        hspec2 = ps.PlottingSpectra(zz, halo10, label="MVELS")
        plot_check(hspec,hspec2,"box", zz)

def test_atten():
    """Plot the effect of the self-shielding correction"""
    halo = myname.get_name(7)
    hspec = ps.PlottingSpectra(3, halo)
    hspec2 = ps.PlottingSpectra(3, halo,savefile="grid_spectra_DLA_no_atten.hdf5")
    plot_check(hspec,hspec2,"no_atten")

def test_noise():
    """Plot the effect of noise on the spectrum"""
    halo = myname.get_name(7)
    hspec = ps.PlottingSpectra(3, halo, snr=0.05)
    hspec2 = ps.PlottingSpectra(3, halo, snr = 0.)
    plot_check(hspec,hspec2,"noise")

    #Higher resolution spectrum
def plot_check(hspec, hspec2, ofile, snap=3):
    """Plot velocity widths for two halos both absolutely and relatively"""
    hspec.plot_vel_width("Si",2, color="red")
    hspec2.plot_vel_width("Si", 2, color="blue")
    vel_data.plot_prochaska_2008_data()
    plt.xlim(1, 1000)
    save_figure(path.join(outdir,"cosmo_vel_width_"+ofile+"_z"+str(snap)))
    plt.clf()
    (vbin,one) = hspec.vel_width_hist("Si",2)
    (vbin,two) = hspec2.vel_width_hist("Si",2)
    if np.size(one) != np.size(two):
        maxx = np.min([np.size(one),np.size(two)])
        plt.semilogx(vbin[:maxx],one[:maxx]/two[:maxx])
    else:
        plt.semilogx(vbin,one/two)
    save_figure(path.join(outdir,"cosmo_rel_vel_width_"+ofile+"_z"+str(snap)))
    plt.clf()
    hspec.plot_f_peak("Si", 2, color="red")
    hspec2.plot_f_peak("Si", 2, color="blue")
    vel_data.plot_extra_stat_hist(True)
    save_figure(path.join(outdir,"cosmo_fpeak_"+ofile))
    plt.clf()
    hspec.plot_Z_vs_vel_width(color="red", color2="darkred")
    hspec2.plot_Z_vs_vel_width(color="blue", color2="purple")
    vel_data.plot_prochaska_2008_correlation(zrange[snap])
    save_figure(path.join(outdir,"cosmo_correlation_"+ofile+"_z"+str(snap)))
    plt.clf()
    hspec.plot_eq_width("Si", 2, 1526, color="red")
    hspec2.plot_eq_width("Si", 2, 1526, color="blue")
    vel_data.plot_si1526_eqw(zrange[snap], nv_table=7)
    save_figure(path.join(outdir,"cosmo_eqwidth_"+ofile+"_z"+str(snap)))
    plt.clf()

def test_tescari_halos(sim, snap):
    """Plot velocity width for spectra through the center of halos, like in Tescari 2009"""
    halo = myname.get_name(sim)
    hspec = ps.PlottingSpectra(snap, halo, savefile="halo_spectra.hdf5")
    hspec.plot_vel_width("Si", 2, color="red")
    vel_data.plot_prochaska_2008_data()
    save_figure(path.join(outdir,"cosmo_tescari_halos"))
    plt.clf()


if __name__ == "__main__":
    test_tescari_halos(7,3)
    test_noise()
    test_atten()
    test_spec_resolution()
    test_box_resolution()
#     test_pecvel()
    test_tophat()

#     plot_vel_width_SiII(0, 3)

#     plot_vel_width_metcol(0,3)
    plot_metal_ion_corr(0,3)
