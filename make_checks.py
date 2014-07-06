"""Make some plots of the velocity widths from the cosmo runs"""
#!/usr/bin env python
# -*- coding: utf-8 -*-

import matplotlib
matplotlib.use('PDF')

import matplotlib.pyplot as plt

import vw_plotspectra as ps
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
    hspec = ps.VWPlotSpectra(snap, halo)
    hspec.plot_metallicity(color="red", ls="-")
    hspec.plot_species_metallicity(species, ion, color="blue", ls="-")
    vel_data.plot_alpha_metal_data((3.5,2.5))
    save_figure(path.join(outdir, "cosmo_metallicity"+str(sim)+"_ion_corr"+str(snap)))
    plt.clf()
    hspec.plot_ion_corr(species, ion)
    save_figure(path.join(outdir, "cosmo_metallicity"+str(sim)+"_rel_ion_corr"+str(snap)))
    plt.clf()

def plot_vel_width_SiII(sim, snap):
    """
       Plot the change in velocity widths between the full calculation and
       setting n(Si+)/n(Si) = n(HI)/n(H)
    """
    halo = myname.get_name(sim)
    #Load from a save file only
    hspec = ps.VWPlotSpectra(snap, halo)
    hspecSi = ps.VWPlotSpectra(snap, halo,savefile="SiHI_spectra.hdf5")
    plot_check(hspec, hspecSi,"SiHI")

def plot_vel_width_SiII_keating(sim, snap):
    """
       Plot the change in velocity widths between the full calculation and
       setting n(Si+)/n(Si) = n(HI)/n(H)
    """
    halo = myname.get_name(sim)
    #Load from a save file only
    hspec = ps.VWPlotSpectra(snap, halo)
    hspecSi2 = ps.VWPlotSpectra(snap, halo,savefile="si_colden_spectra.hdf5")
    plot_check(hspec, hspecSi2,"SiHI_keating")

def test_spec_resolution():
    """Plot the velocity widths for different spectral resolutions"""
    halo = myname.get_name(7)
    #Higher resolution spectrum
    hspec = ps.VWPlotSpectra(3, halo, savefile="grid_spectra_DLA.hdf5")
    hspec2 = ps.VWPlotSpectra(3, halo, savefile="grid_spectra_DLA_res.hdf5")
    plot_check(hspec,hspec2,"specres")

def test_vel_abswidth():
    """Plot the velocity widths for different minimum absorber widths"""
    halo = myname.get_name(7)
    #Higher resolution spectrum
    hspec = ps.VWPlotSpectra(3, halo)
    hspec2 = ps.VWPlotSpectra(3, halo)
    hspec2.minwidth = 250.
    plot_check(hspec,hspec2,"abswidth")

def test_pecvel():
    """Plot the velocity widths with and without peculiar velocities"""
    halo = myname.get_name(7)
    #Higher resolution spectrum
    hspec = ps.VWPlotSpectra(3, halo, savefile="grid_spectra_DLA.hdf5")
    hspec2 = ps.VWPlotSpectra(3, halo, None, None, savefile="grid_spectra_DLA_pecvel.hdf5")
    plot_check(hspec,hspec2,"pecvel")

def test_tophat():
    """Plot the velocity widths with and with top hat vs SPH"""
    halo = myname.get_name(7)
    #Higher resolution spectrum
    hspec = ps.VWPlotSpectra(3, halo, savefile="grid_spectra_DLA.hdf5")
    hspec2 = ps.VWPlotSpectra(3, halo, None, None, savefile="grid_spectra_DLA_tophat.hdf5")
    plot_check(hspec,hspec2,"tophat")

def test_lowres():
    """Plot the velocity widths with and with top hat vs SPH"""
    halo = myname.get_name(0)
    halolow = myname.get_name(0, ff=False)
    #Higher resolution spectrum
    hspec = ps.VWPlotSpectra(3, halo, savefile="grid_spectra_DLA.hdf5")
    hspec2 = ps.VWPlotSpectra(60, halolow, None, None, savefile="rand_spectra_DLA.hdf5")
    plot_check(hspec,hspec2,"lowres")

def test_box_resolution():
    """Plot the velocity widths for different size boxes"""
    halo = myname.get_name(7)
    halo10 = myname.get_name(7,box=7.5)
#     for zz in (1,3,5):
    zz = 3
    hspec = ps.VWPlotSpectra(zz, halo, label="DEF")
    hspec2 = ps.VWPlotSpectra(zz, halo10, label="SMALL")
    plot_check(hspec,hspec2,"box", zz)

def test_big_box():
    """Plot the velocity widths for different size boxes"""
    halo = myname.get_name(0)
    halobig = path.expanduser("~/data/Illustris")
    hspec = ps.VWPlotSpectra(3, halo, label="DEF")
    hspec2 = ps.VWPlotSpectra(59, halobig, label="ILLUS")
    plot_check(hspec,hspec2,"bigbox")

def test_gfm_shield():
    """Plot the velocity widths for dynamical self-shielding vs post-processed self-shielding."""
    halo = myname.get_name(7)
    halo2 = myname.get_name('B')
    hspec = ps.VWPlotSpectra(3, halo, label="2xUV")
    hspec2 = ps.VWPlotSpectra(3, halo2, label="NOSHIELD")
    plot_check(hspec,hspec2,"gfm_shield")
    hspec = ps.VWPlotSpectra(5, halo, label="2xUV")
    hspec2 = ps.VWPlotSpectra(5, halo2, label="NOSHIELD")
    plot_check(hspec,hspec2,"gfm_shield", snap=5)

class NoFilt(ps.VWPlotSpectra):
    def get_filt(self, elem, ion):
        return ps.VWPlotSpectra.get_filt(self, elem, ion, 100)

def test_filt():
    """Plot impact of filtering low-metallicity systems."""
    halo = myname.get_name(7)
    hspec = ps.VWPlotSpectra(3, halo, label="FILT")
    hspec2 = NoFilt(3, halo, label="NOFILT")
    plot_check(hspec,hspec2,"filtering")

def test_atten():
    """Plot the effect of the self-shielding correction"""
    halo = myname.get_name(7)
    hspec = ps.VWPlotSpectra(3, halo, label="ATTEN")
    hspec2 = ps.VWPlotSpectra(3, halo,savefile="grid_spectra_DLA_no_atten.hdf5",label="NOATTEN")
    plot_check(hspec,hspec2,"no_atten")

def test_shield():
    """Plot velocity width for spectra using self-shielding like in Tescari 2009"""
    halo = myname.get_name(7)
    hspec = ps.VWPlotSpectra(3, halo)
    hspec2 = ps.VWPlotSpectra(3, halo,savefile="grid_spectra_DLA_noshield.hdf5")
    plot_check(hspec,hspec2,"no_shield")

def test_noise():
    """Plot the effect of noise on the spectrum"""
    halo = myname.get_name(7)
    hspec = ps.VWPlotSpectra(3, halo, snr=0.,label="No Noise")
    hspec2 = ps.VWPlotSpectra(3, halo, snr = 20.,label="Noise")
    plot_check(hspec,hspec2,"noise")

    #Higher resolution spectrum
def plot_check(hspec, hspec2, ofile, snap=3):
    """Plot velocity widths for two halos both absolutely and relatively"""
    hspec.plot_vel_width("Si",2, color="red")
    hspec2.plot_vel_width("Si", 2, color="blue", ls="--")
    plt.legend()
    vel_data.plot_prochaska_2008_data()
    plt.xlabel(r"$v_\mathrm{90}$ (km s$^{-1}$)")
    plt.xlim(10, 1000)
    save_figure(path.join(outdir,"cosmo_vel_width_"+ofile+"_z"+str(snap)))
    plt.clf()
    (vbin,one) = hspec.vel_width_hist("Si",2)
    (vbin,two) = hspec2.vel_width_hist("Si",2)
    if np.size(one) != np.size(two):
        maxx = np.min([np.size(one),np.size(two)])
        plt.semilogx(vbin[:maxx],one[:maxx]/two[:maxx])
    else:
        plt.semilogx(vbin,one/two)
    plt.legend()
    save_figure(path.join(outdir,"cosmo_rel_vel_width_"+ofile+"_z"+str(snap)))
    plt.clf()
    hspec.plot_f_peak("Si", 2, color="red")
    hspec2.plot_f_peak("Si", 2, color="blue", ls="--")
    vel_data.plot_extra_stat_hist(True)
    plt.legend()
    save_figure(path.join(outdir,"cosmo_fpeak_"+ofile))
    plt.clf()
    hspec.plot_Z_vs_vel_width(color="red", color2="darkred")
    hspec2.plot_Z_vs_vel_width(color="blue", color2="purple")
    vel_data.plot_prochaska_2008_correlation(zrange[snap])
    plt.legend()
    save_figure(path.join(outdir,"cosmo_correlation_"+ofile+"_z"+str(snap)))
    plt.clf()
    hspec.plot_eq_width("Si", 2, 1526, color="red")
    hspec2.plot_eq_width("Si", 2, 1526, color="blue", ls="--")
    vel_data.plot_si1526_eqw(zrange[snap], nv_table=7)
    plt.legend()
    save_figure(path.join(outdir,"cosmo_eqwidth_"+ofile+"_z"+str(snap)))
    plt.clf()

def test_tescari_halos(sim, snap):
    """Plot velocity width for spectra through the center of halos, like in Tescari 2009"""
    halo = myname.get_name(sim, box=10)
    hspec = ps.VWPlotSpectra(snap, halo, savefile="halo_spectra_2.hdf5",cdir=path.expanduser("~/codes/cloudy_tables/ion_out_no_atten/"))
    hspec.plot_vel_width("Si", 2, color="red")
    vel_data.plot_prochaska_2008_data()
    save_figure(path.join(outdir,"cosmo_tescari_halos"))
    plt.clf()

if __name__ == "__main__":
#     test_shield()
    test_vel_abswidth()
    test_gfm_shield()
    test_tescari_halos(5,3)
    test_noise()
    test_atten()
    test_spec_resolution()
    test_lowres()
    test_big_box()
    test_box_resolution()
#     test_pecvel()
    test_tophat()

    plot_vel_width_SiII(7, 3)
    plot_vel_width_SiII_keating(7, 3)

#     plot_vel_width_metcol(0,3)
    plot_metal_ion_corr(7,3)
    plot_metal_ion_corr(0,3)
