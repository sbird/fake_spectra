# -*- coding: utf-8 -*-
"""Contains the plotting-specific functions for the spectrum analysis code."""

from __future__ import print_function
import numpy as np
#import leastsq as lsq
import matplotlib.pyplot as plt

from . import spectra

try:
    xrange(1)
except NameError:
    xrange = range

class PlottingSpectra(spectra.Spectra):
    """Class to plot things connected with spectra."""
    def __init__(self,num, base, cofm=None, axis=None, label="", snr=0., load_halo=True,**kwargs):
        spectra.Spectra.__init__(self,num, base, cofm=cofm, axis=axis, snr=snr, load_halo=load_halo, **kwargs)
        self.label=label

    def plot_eq_width(self, elem, ion, line, dv=0.1, color="red", ls="-"):
        """Plot the equivalent width histogram of this snapshot
        Parameters:
            elem - element to use
            ion - ionisation state: 1 is neutral.
            line - line number to use
            dv - bin spacing
        """
        (vbin, eqw) = self.eq_width_hist(elem, ion, line, dv)
        plt.plot(vbin, eqw, color=color, lw=3, ls=ls,label=self.label)

    def plot_spectrum(self, elem, ion, line, spec_num, flux=True, xlims=(-500,500), color="blue",ls="-", offset=0):
        """Plot an spectrum, centered on the maximum tau.
        Parameters:
            elem, ion, line - line profile to plot.
            spec_num - number of the spectrum, in the catalogue, to plot
            flux - if False, print tau, if True print e^-tau
        """
        tau = self.get_tau(elem, ion, line, spec_num, noise=True)
        peak = np.where(tau == np.max(tau))[0][0]
        szt = int(np.size(tau)/2)
        tau_l = np.roll(tau, szt - peak+int(offset*self.dvbin))
        xaxis = (np.arange(0,np.size(tau))-szt)*self.dvbin
        self.plot_spectrum_raw(tau_l,xaxis, xlims, flux, color=color,ls=ls)
        return peak

    def plot_spectrum_raw(self, tau,xaxis,xlims, flux=True, color="blue",ls="-"):
        """Plot an array of optical depths, centered on the largest point,
           and marking the 90% velocity width.
           offset: offset in km/s for the x-axis labels"""
        #Make sure we were handed a single spectrum
        assert np.size(np.shape(tau)) == 1
        if flux:
            plt.plot(xaxis,np.exp(-tau), color=color,ls=ls)
        else:
            plt.plot(xaxis,tau,color=color, ls=ls)
        plt.xlim(xlims)
        plt.xlabel(r"v (km s$^{-1}$)")
        if flux:
            plt.ylabel(r"$\mathcal{F}$")
            plt.ylim(-0.05,1.05)
        else:
            plt.ylabel(r"$\tau$")
            plt.ylim(-0.1,np.max(tau)+0.2)
        return xaxis[0]

    def plot_density(self, elem, ion, num, thresh=1e-9, color="blue"):
        """Plot the density of an ion along a sightline"""
        den = self.get_density(elem, ion)
        mcol = np.max(den[num])
        ind_m = np.where(den[num] == mcol)[0][0]
        den = np.roll(den[num], int(np.size(den[num])/2) - ind_m)
        phys = self.dvbin/self.velfac
        #Add one to avoid zeros on the log plot
        plt.semilogy(np.arange(0,np.size(den))*phys-np.size(den)/2*phys,den+1e-30, color=color)
        plt.xlabel(r"x (kpc h$^{-1}$)")
        plt.ylabel(r"n (cm$^{-3}$)")
        #Set limits
        ind = np.where(den > thresh)
        if np.size(ind) > 0:
            dxlim = np.max(np.abs((ind[0][0]*phys-np.size(den)/2*phys-50, ind[0][-1]*phys-np.size(den)/2*phys+50)))
            plt.xlim(-1.*dxlim,dxlim)
        else:
            dxlim = np.size(den)*phys
        return dxlim

    def plot_den_to_tau(self, elem, ion, num, thresh = 1e-10,xlim=100, voff = 0., xscale=1):
        """Make a plot connecting density on the low x axis to optical depth on the high x axis.
        Arguments:
            elem, ion - ionic species to plot
            num - index of spectrum shown
            thresh - density threshold above with to track the pixels
            xlim - width of shown plot in km/s
            voff - constant value to shift the high x axis by."""
        #Number of points to draw for each line
        npix = 10
        #Get densities above threshold
        den = self.get_density(elem, ion)[num]
        imax = np.where(den == np.max(den))[0][0]
        ind = np.where(den > thresh)
        #Get peculiar velocity along sightline
        ax = self.axis[num]-1
        vel = self.get_velocity(elem, ion)[num, :, ax]
        #Adjust the axis offset.
        vel -= vel[imax]-voff
        #Convert pixel coordinates to offsets from peak
        ind = np.ravel(ind)
        coord = (ind - imax)*self.dvbin
        coord[np.where(coord > self.nbins/2)] = coord - self.nbins
        coord[np.where(coord < -self.nbins/2)] = coord + self.nbins
        for (cc, ii) in zip(coord, ind):
            x = np.linspace(cc/xscale,cc+vel[ii], npix)
            y = np.linspace(0,1,npix)
            plt.plot(x,y,ls="-", color="black")
        plt.xlabel(r"v (km s$^{-1}$)")
        plt.xlim(-1.*xlim, xlim)
        plt.ylim(0,1)
        plt.yticks(())

    def plot_temp(self, elem, ion):
        """Make a contour plot for the density weighted temperature for each spectrum"""
        temp = self.get_temp(elem, ion)
        den = self.get_density(elem, ion)
        ind = np.where(den < 1e-6)
        den[ind] = 0.
        temps = np.sum(temp*den, axis=1)/np.sum(den, axis=1)
        ind2 = np.where(temps < 1e5)
        print(np.median(temps)," filt: ", np.median(temps[ind2]))
        self._plot_2d_contour(np.sum(den,axis=1)[ind2], temps[ind2], 40, name="Temp Density", color="blue", color2="darkblue", ylog=False, xlog=True, fit=False, sample = 300)
        plt.xlabel(r"n (cm$^{-3}$)")
        plt.ylabel(r"T (K)")
        plt.ylim(0,2e4)

    def plot_cddf(self,elem = "H", ion = 1, dlogN=0.2, minN=13, maxN=23., color="blue", moment=False, dX=True):
        """Plots the column density distribution function. """
        (NHI,f_N)=self.column_density_function(elem, ion, dlogN,minN-1,maxN+1,dX=dX)
        if moment:
            f_N *= NHI
        plt.loglog(NHI,f_N,color=color, label=self.label)
        ax=plt.gca()
        ax.set_xlabel(r"$N (\mathrm{cm}^{-2})$")
        ax.set_ylabel(r"$f(N) (\mathrm{cm}^2)$")
        plt.xlim(10**minN, 10**maxN)
#         if moment:
#             plt.ylim(1e-4,1)

    def _plot_metallicity(self, met, nbins=20,color="blue", ls="-"):
        """Plot the distribution of metallicities"""
        bins=np.linspace(-3,0,nbins)
        mbin = np.array([(bins[i]+bins[i+1])/2. for i in range(0,np.size(bins)-1)])
        #Abs. distance for entire spectrum
        hist = np.histogram(np.log10(met),bins,density=True)[0]
        plt.plot(mbin,hist,color=color,label=self.label,ls=ls)

    def plot_metallicity(self, nbins=20,color="blue", ls="-", width=0.):
        """Plot the distribution of metallicities"""
        met = self.get_metallicity(width=width)
        self._plot_metallicity(met,nbins,color,ls)

    def plot_species_metallicity(self, species, ion, nbins=20,color="blue", ls="-"):
        """Plot the distribution of metallicities from an ionic species"""
        met = self.get_ion_metallicity(species,ion)
        self._plot_metallicity(met,nbins,color,ls)

    def plot_ion_corr(self, species, ion, nbins=80,color="blue",ls="-"):
        """Plot the difference between the single-species ionisation and the metallicity from GFM_Metallicity"""
        met = np.log10(self.get_metallicity())
        ion_met = np.log10(self.get_ion_metallicity(species, ion))
        diff = 10**(ion_met - met)
        print(np.max(diff), np.min(diff), np.median(diff))
        bins=np.linspace(-1,1,nbins)
        mbin = np.array([(bins[i]+bins[i+1])/2. for i in range(0,np.size(bins)-1)])
        hist = np.histogram(np.log10(diff),bins,density=True)[0]
        plt.plot(mbin,hist,color=color,label=self.label,ls=ls)

    def plot_eq_width_vs_col_den(self, elem, ion, line):
        """Plot the equivalent width vs the column density along the sightline for each spectrum."""
        eqw = np.log10(self.equivalent_width(elem, ion, line))
        colden = np.sum(self.get_col_density(elem, ion), axis=1)
        plt.semilogy(eqw, colden,'o')
        plt.xlabel(r"W $(\AA)$")
        plt.ylabel(r"N$_\mathrm{HI}$ (cm$^{2}$)")

    def _plot_2d_contour(self, xvals, yvals, nbins, name="x y", color="blue", color2="darkblue", ylog=True, xlog=True, fit=False, sample=40.):
        """Helper function to make a 2D contour map of a correlation, as well as the best-fit linear fit"""
        if ylog:
            yvals = np.log10(yvals)
        if xlog:
            xvals = np.log10(xvals)
        (H, xedges, yedges) = np.histogram2d(xvals, yvals,bins=nbins)
        xbins=np.array([(xedges[i+1]+xedges[i])/2 for i in xrange(0,np.size(xedges)-1)])
        ybins=np.array([(yedges[i+1]+yedges[i])/2 for i in xrange(0,np.size(yedges)-1)])
        xx = np.logspace(np.min(xbins), np.max(xbins),15)
        ax = plt.gca()
        if ylog:
            ybins = 10**ybins
            ax.set_yscale('log')
        if xlog:
            xbins = 10**xbins
            ax.set_xscale('log')
        plt.contourf(xbins,ybins,H.T,(self.NumLos/sample)*np.array([0.15,1,10]),colors=(color,color2,"black"),alpha=0.5)
        #if fit:
        #    (intercept, slope, _) = lsq.leastsq(xvals,yvals)
        #    plt.loglog(xx, 10**intercept*xx**slope, color="black",label=self.label, ls="--")
