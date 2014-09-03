# -*- coding: utf-8 -*-
"""Contains the plotting-specific functions for the spectrum analysis code."""

import spectra
import numpy as np
import leastsq as lsq
import matplotlib.pyplot as plt

class PlottingSpectra(spectra.Spectra):
    """Class to plot things connected with spectra."""
    def __init__(self,num, base, cofm=None, axis=None, res=1., savefile="grid_spectra_DLA.hdf5",label="", snr=0., spec_res = 8., cdir=None):
        spectra.Spectra.__init__(self,num, base, cofm, axis, res, savefile=savefile, snr=snr, spec_res=spec_res, cdir=cdir)
        self.label=label

    def plot_eq_width(self, elem, ion, line, dv=0.1, eq_cut = 0.002, color="red", ls="-"):
        """Plot the velocity widths of this snapshot
        Parameters:
            elem - element to use
            ion - ionisation state: 1 is neutral.
            line - line number to use
            dv - bin spacing
        """
        (vbin, eqw) = self.eq_width_hist(elem, ion, line, dv, eq_cut=eq_cut)
        plt.plot(vbin, eqw, color=color, lw=3, ls=ls,label=self.label)

    def plot_spectrum(self, elem, ion, line, num, flux=True, xlims=(-500,500), color="blue"):
        """Plot an spectrum, centered on the maximum tau"""
        tau = self.get_tau(elem, ion, line, num, noise=True)
        peak = np.where(tau == np.max(tau))[0]
        tau_l = np.roll(tau, peak)
        xaxis = np.arange(-np.size(tau)/2,np.size(tau)/2)*self.dvbin
        self.plot_spectrum_raw(tau_l,xaxis, xlims, flux, color=color)
        return peak

    def plot_spectrum_raw(self, tau,xaxis,xlims, flux=True, color="blue"):
        """Plot an array of optical depths, centered on the largest point,
           and marking the 90% velocity width.
           offset: offset in km/s for the x-axis labels"""
        #Make sure we were handed a single spectrum
        assert np.size(np.shape(tau)) == 1
        if flux:
            plt.plot(xaxis,np.exp(-tau), color=color)
        else:
            plt.plot(xaxis,tau,color=color)
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
        den = np.roll(den[num], np.size(den[num])/2 - ind_m)
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
        print np.median(temps)," filt: ", np.median(temps[ind2])
        self._plot_2d_contour(np.sum(den,axis=1)[ind2], temps[ind2], 40, name="Temp Density", color="blue", color2="darkblue", ylog=False, xlog=True, fit=False, sample = 300)
        plt.xlabel(r"n (cm$^{-3}$)")
        plt.ylabel(r"T (K)")
        plt.ylim(0,2e4)

    def plot_cddf(self,elem = "H", ion = 1, dlogN=0.2, minN=13, maxN=23., color="blue", moment=False):
        """Plots the column density distribution function. """
        (NHI,f_N)=self.column_density_function(elem, ion, dlogN,minN-1,maxN+1)
        if moment:
            f_N *= NHI
        plt.loglog(NHI,f_N,color=color, lw = 3)
        ax=plt.gca()
        ax.set_xlabel(r"$N_\mathrm{HI} (\mathrm{cm}^{-2})$")
        ax.set_ylabel(r"$f(N) (\mathrm{cm}^2)$")
        plt.xlim(10**minN, 10**maxN)
        if moment:
            plt.ylim(1e-4,1)

    def _plot_breakdown(self, array, filt, low, high, labels, dv, log=True):
        """
        Helper function to plot something broken down by halo mass
        """
        #Find virial velocity
        (halo, _) = self.find_nearest_halo()
        ind = np.where(halo[filt] > 0)
        virial = self.virial_vel(halo[filt][ind])
        array = array[filt]
        #Make bins
        if log:
            func = plt.semilogx
            v_table = 10**np.arange(np.min(np.log10(array)),np.max(np.log10(array)) , dv)
        else:
            func = plt.plot
            v_table = np.arange(np.min(array),np.max(array) , dv)
        vbin = np.array([(v_table[i]+v_table[i+1])/2. for i in range(0,np.size(v_table)-1)])
        #Histogram of vel width
        vhist = np.histogram(array, v_table)[0]
        vhist[np.where(vhist == 0)] = 1
        colors = ("red", "purple", "cyan")
        lss = ("--", ":", "-")
        #Histogram of vel width for all halos in given virial velocity bin
        for ii in xrange(len(low)):
            vind = np.where((virial > low[ii])*(virial < high[ii]))
            vhist2 = np.histogram(array[ind][vind], v_table)[0]
            func(vbin, vhist2/(1.*vhist), color=colors[ii], ls=lss[ii], label=labels[ii])
#         vind = np.where(halo[filt] < 0)
#         vhist2 = np.histogram(array[vind], v_table)[0]
#         func(vbin, vhist2/(1.*vhist), color="grey", ls="-.", label="Field")

    def _plot_metallicity(self, met, nbins=20,color="blue", ls="-"):
        """Plot the distribution of metallicities"""
        bins=np.linspace(-3,0,nbins)
        mbin = np.array([(bins[i]+bins[i+1])/2. for i in range(0,np.size(bins)-1)])
        #Abs. distance for entire spectrum
        hist = np.histogram(np.log10(met),bins,density=True)[0]
        plt.plot(mbin,hist,color=color,label=self.label,ls=ls)

    def plot_metallicity(self, nbins=20,color="blue", ls="-"):
        """Plot the distribution of metallicities"""
        met = self.get_metallicity()
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
        print np.max(diff), np.min(diff), np.median(diff)
        bins=np.linspace(-1,1,nbins)
        mbin = np.array([(bins[i]+bins[i+1])/2. for i in range(0,np.size(bins)-1)])
        hist = np.histogram(np.log10(diff),bins,density=True)[0]
        plt.plot(mbin,hist,color=color,label=self.label,ls=ls)

    def plot_Z_vs_mass(self,color="blue", color2="darkblue"):
        """Plot the correlation between mass and metallicity, with a fit"""
        (halo, _) = self.find_nearest_halo()
        ind = np.where(halo > 0)
        met = self.get_metallicity()[ind]
        mind = np.where(met > 1e-4)
        halo = halo[ind]
        mass = self.sub_mass[halo]
        mass = mass[mind]
        met = met[mind]
        self._plot_2d_contour(mass+0.1, met, 10, "Z mass", color, color2)
        plt.ylim(1e-4,1)

    def _plot_xx_vs_mass(self, xx, name = "xx", color="blue", color2="darkblue", log=True):
        """Helper function to plot something against virial velocity"""
        (halo, _) = self.find_nearest_halo()
        ii = self.get_filt("Si",2)
        ind = np.where(halo[ii] > 0)
        halo = halo[ii][ind]
        xx = xx[ii][ind]
        virial = self.virial_vel(halo)+0.1
        self._plot_2d_contour(virial, xx, 10, name+" virial velocity", color, color2, ylog=log)

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
        if fit:
            (intercept, slope, _) = lsq.leastsq(xvals,yvals)
            plt.loglog(xx, 10**intercept*xx**slope, color="black",label=self.label, ls="--")
