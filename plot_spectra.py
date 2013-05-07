# -*- coding: utf-8 -*-
"""Contains the plotting-specific functions for the spectrum analysis code."""

import convert_cloudy
import spectra
import halospectra
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt

class PlottingSpectra(spectra.Spectra):
    """Class to plot things connected with spectra."""
    def __init__(self,num, base, cofm, axis, res=1., savefile="rand_spectra_DLA.hdf5"):
        spectra.Spectra.__init__(self,num, base, cofm, axis, res, savefile)

    def plot_vel_width(self, elem, line, dv=0.1, HI_cut = None, met_cut = 1e13, unres = 10, color="red"):
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
        plt.loglog(vbin, vels, color=color, lw=3)

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
        (low, high) = self._vel_width_bound(tau_l, tot_tau)
        if high - low > 0:
            plt.plot([low,low],[0,1])
            plt.plot([high,high],[0,1])
        plt.text(high+self.dvbin*30,0.5,r"$\delta v_{90} = "+str(np.round(high-low,1))+r"$")
        plt.ylim(-0.05,1.05)
        plt.xlim(0,np.size(tau_l)*self.dvbin)
        return (low, high)

    def plot_spectrum_density_velocity(self, elem, ion, i):
        """Plot the spectrum of a line, centered on the deepest point,
           and marking the 90% velocity width."""
        #  Size of a single velocity bin
        tau = self.get_tau(elem, ion)
        col_den = self.get_col_density(elem, ion)[i,:]
        #Deal with periodicity by making sure the deepest point is in the middle
        tau_l = tau[i,:]
        tmax = np.max(tau_l)
        ind_m = np.where(tau_l == tmax)[0][0]
        tau_l = np.roll(tau_l, np.size(tau_l)/2- ind_m)
        col_den = np.roll(col_den, np.size(tau_l)/2- ind_m)
        plt.subplot(311)
        (low,high) = self.plot_spectrum(tau, i)
        plt.xlim(low-50,high+50)
        plt.subplot(312)
        ind = np.where(col_den > 1)
        plt.semilogy(np.arange(0,np.size(tau_l))[ind]*self.dvbin,col_den[ind])
        #plt.xlim(0,np.size(tau_l)*self.dvbin)
        plt.xlim(low-50,high+50)
        plt.ylim(1e5,1e20)
        plt.yticks(np.array([1e10, 1e15,1e20]))
        plt.subplot(313)
        vel = self.get_vel(elem, ion)[i,:]
        vel = np.roll(vel, np.size(tau_l)/2- ind_m)
        ind = np.where(np.abs(vel) > 1)
        plt.semilogy(np.arange(0,np.size(tau_l))[ind]*self.dvbin,np.abs(vel[ind]))
        #plt.xlim(0,np.size(tau_l)*self.dvbin)
        plt.xlim(low-50,high+50)

    def plot_col_density(self, elem, ion):
        """Plot the maximal column density in each sightline against vel_width, assuming rho and tau were already calculated"""
        col_dens = self.get_col_density(elem, ion)
        vels = self.vel_width(self.metals[(elem, ion)][3])
        plt.loglog(np.max(col_dens,axis=1),vels)

    def plot_cddf(self,elem = "H", ion = 1, dlogN=0.2, minN=13, maxN=23., color="blue"):
        """Plots the column density distribution function. """
        (NHI,f_N)=self.column_density_function(elem, ion, dlogN,minN-1,maxN+1)
        plt.loglog(NHI,f_N,color=color, lw = 3)
        ax=plt.gca()
        ax.set_xlabel(r"$N_\mathrm{HI} (\mathrm{cm}^{-2})$")
        ax.set_ylabel(r"$f(N) (\mathrm{cm}^2)$")
        plt.xlim(10**minN, 10**maxN)
        plt.ylim(1e-26,1e-18)

    def plot_sep_frac(self,elem = "Si", ion = 2, thresh = 1e-2, mindist = 15, dv = 0.1):
        """Plots the fraction of spectra in each velocity width bin which are separated.
        Threshold is as a percentage of the maximum value.
        mindist is in km/s
        """
        sep = self.get_separated(elem, ion, thresh,mindist)
        vels = self.vel_width(self.get_tau(elem, ion))
        v_table = 10**np.arange(0, np.log10(np.max(vels)), dv)
        vbin = np.array([(v_table[i]+v_table[i+1])/2. for i in range(0,np.size(v_table)-1)])
        hist1 = np.histogram(vels, v_table)
        hist2 = np.histogram(vels[sep],v_table)
        hist1[0][np.where(hist1[0] == 0)] = 1
        plt.semilogx(vbin, hist2[0]/(1.*hist1[0]))

    def plot_metallicity(self, nbins=20,color="blue"):
        """Plot the distribution of metallicities"""
        bins=np.logspace(-3,0,nbins)
        mbin = np.array([(bins[i]+bins[i+1])/2. for i in range(0,np.size(bins)-1)])
        met = self.get_metallicity()
        hist = np.histogram(met,bins)[0]
        #Abs. distance for entire spectrum
        hist = np.histogram(np.log10(met),np.log10(bins),density=True)[0]
        plt.semilogx(mbin,hist,color=color)

    def plot_Z_vs_vel_width(self,elem="Si", line=2):
        """Plot the correlation between metallicity and velocity width"""
        met = self.get_metallicity()
        tau = self.get_tau(elem, line, 2)
        ind = self.get_filt(elem, line, 10**20.3, 10**12)
        met = met[ind]
        vel_width = self.vel_width(tau[ind])
        ind2 = np.where(vel_width > 15)
        plt.loglog(vel_width[ind2],met[ind2], 'x')
        (slope, intercept, rval, pval, sd) = stats.linregress(np.log10(vel_width[ind2]),np.log10(met[ind2]))
        print "corr: ",rval, pval, sd
        xx = np.logspace(np.log10(np.min(vel_width)), np.log10(np.max(vel_width)),15)
        plt.loglog(xx,10**intercept*xx**slope)
        plt.xlim(10,2e3)

class PlotHaloSpectra(halospectra.HaloSpectra, PlottingSpectra):
    """Class to plot things connected with spectra."""
    def __init__(self,num, base, repeat = 3, minpart = 400, res = 1., savefile="halo_spectra_DLA.hdf5"):
        halospectra.HaloSpectra.__init__(self,num, base, repeat, minpart, res, savefile)

class PlotIonDensity:
    """Class to plot the ionisation fraction of elements as a function of density"""
    def __init__(self, red):
        self.cloudy_table = convert_cloudy.CloudyTable(red)
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

