# -*- coding: utf-8 -*-
"""Contains the plotting-specific functions for the spectrum analysis code."""

import convert_cloudy
import spectra
import numpy as np
import matplotlib.pyplot as plt

class PlottingSpectra(spectra.Spectra):
    """Class to plot things connected with spectra."""
    def __init__(self,num, base, cofm, axis, nbins = 1024, cloudy_dir="/home/spb/codes/ArepoCoolingTables/tmp_spb/", savefile=None):
        spectra.Spectra.__init__(self,num, base, cofm, axis, nbins, cloudy_dir, savefile)

    def plot_vel_width(self, elem, line, dv=0.1, HI_cut = None, met_cut = 1e13, unres = 20, color="red"):
        """Plot the velocity widths of this snapshot
        Parameters:
            elem - element to use
            line - line to use (the components of this line must be pre-computed and stored in self.metals)
            dv - bin spacing
            HI_cut - Prochaska used a subsample of spectra containing a DLA.
                     If this value is not None, consider only HI column densities above this threshold.
                     If the spectra are taken within the halo virial radius, this does not make much of a difference.
            met_cut - Discard spectra whose maximal metal column density is below this level.
                      Removes unobservable systems.
            unres - Remove systems with velocity widths below this value, where they are affected
                    by the pixel size of the spectra.

        """
        (vbin, vels) = self.vel_width_hist(elem, line, dv, HI_cut, met_cut, unres)
        plt.semilogy(vbin, vels, color=color)

    def plot_spectrum(self, tau, i):
        """Plot the spectrum of a line, centered on the deepest point,
           and marking the 90% velocity width."""
        #  Size of a single velocity bin
        tot_tau = np.sum(tau[i,:])
        #Deal with periodicity by making sure the deepest point is in the middle
        tau_l = tau[i,:]
        tmax = np.max(tau_l)
        ind_m = np.where(tau_l == tmax)[0][0]
        tau_l = np.roll(tau_l, np.size(tau_l)/2- ind_m)
        plt.plot(np.arange(0,np.size(tau_l))*self.dvbin,np.exp(-tau_l))
        cum_tau = np.cumsum(tau_l)
        ind_low = np.where(cum_tau > 0.05 * tot_tau)
        low = ind_low[0][0]*self.dvbin
        ind_high = np.where(cum_tau > 0.95 * tot_tau)
        high = ind_high[0][0]*self.dvbin
        if high - low > 0:
            plt.plot([low,low],[0,1])
            plt.plot([high,high],[0,1])
        plt.text(high+self.dvbin*30,0.5,r"$\delta v_{90} = "+str(np.round(high-low,1))+r"$")
        plt.ylim(-0.05,1.05)
        plt.xlim(0,np.size(tau_l)*self.dvbin)

    def plot_col_density(self, elem, ion):
        """Plot the maximal column density in each sightline against vel_width, assuming rho and tau were already calculated"""
        col_dens = self.get_col_density(elem, ion)
        vels = self.vel_width(self.metals[(elem, ion)][3])
        plt.loglog(np.max(col_dens,axis=1),vels)


class PlotIonDensity:
    """Class to plot the ionisation fraction of elements as a function of density"""
    def __init__(self, red, cloudy_dir="/home/spb/codes/ArepoCoolingTables/tmp_spb/"):
        self.cloudy_table = convert_cloudy.CloudyTable(cloudy_dir, red)
        self.red = red

    def iondensity(self,elem,ion, metal = 0.1, den=(-2.,3)):
        """Plot the ionisation fraction of an ionic species as a function of hydrogen density.
        Arguments:
             elem, ion - specify the species to plot
             metal - metallicity as a fraction of solar for this species
             den - range of densities to plot
        """
        #Bins in density
        dens = 10**np.arange(den[0],den[1],0.2)
        mass_frac = self.cloudy_table.get_solar(elem)*metal*np.ones(np.size(dens))
        ionfrac = self.cloudy_table.ion(elem, ion, mass_frac, dens)
        plt.loglog(dens,ionfrac)

