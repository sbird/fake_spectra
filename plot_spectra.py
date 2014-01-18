# -*- coding: utf-8 -*-
"""Contains the plotting-specific functions for the spectrum analysis code."""

import convert_cloudy
import spectra
import halospectra
import numpy as np
import leastsq as lsq
import kstest as ks
import matplotlib.pyplot as plt

class PlottingSpectra(spectra.Spectra):
    """Class to plot things connected with spectra."""
    def __init__(self,num, base, cofm=None, axis=None, res=1., savefile="grid_spectra_DLA.hdf5",label=""):
        spectra.Spectra.__init__(self,num, base, cofm, axis, res, savefile=savefile)
        self.label=label

    def plot_vel_width(self, elem, ion, dv=0.1, met_cut = 1e13, color="red", ls="-"):
        """Plot the velocity widths of this snapshot
        Parameters:
            elem - element to use
            ion - ionisation state: 1 is neutral.
            dv - bin spacing
            met_cut - Discard spectra whose maximal metal column density is below this level.
                      Removes unobservable systems.

        """
        (vbin, vels) = self.vel_width_hist(elem, ion, dv, met_cut=met_cut)
        plt.semilogx(vbin, vels, color=color, lw=3, ls=ls,label=self.label)

    def plot_f_meanmedian(self, elem, ion, dv=0.03, met_cut = 1e13, color="red", ls="-"):
        """
        Plot an f_mean_median histogram
        For args see plot_vel_width
        """
        (vbin, vels) = self.f_meanmedian_hist(elem, ion, dv, met_cut=met_cut)
        plt.plot(vbin, vels, color=color, lw=3, ls=ls,label=self.label)

    def plot_f_peak(self, elem, ion, dv=0.03, met_cut = 1e13, color="red", ls="-"):
        """
        Plot an f_peak histogram
        For args see plot_vel_width
        """
        (vbin, vels) = self.f_peak_hist(elem, ion, dv, met_cut=met_cut)
        plt.plot(vbin, vels, color=color, lw=3, ls=ls,label=self.label)

    def plot_equivalent_width(self, elem, ion, line, dv=0.1, color="red", ls="-"):
        """Plot the equivalent widths a spectrum.
        Parameters:
            elem - element to use
            ion - ionisation state: 1 is neutral.
            line - transition to plot equivalent width
            dv - bin spacing
        """
        ww = self.equivalent_width(elem, ion, line)
        w_table = 10**np.arange(np.log10(np.min(ww)), np.log10(np.max(ww)), dv)
        wbin = np.array([(w_table[i]+w_table[i+1])/2. for i in range(0,np.size(w_table)-1)])
        whist = np.histogram(np.log10(ww),np.log10(w_table), density=True)[0]
        plt.semilogx(wbin, whist, color=color, lw=3, ls=ls)

    def plot_spectrum(self, tau, width=None):
        """Plot the spectrum of a line, centered on the deepest point,
           and marking the 90% velocity width."""
        #  Size of a single velocity bin
        tot_tau = np.sum(tau)
        #Deal with periodicity by making sure the deepest point is in the middle
        tau_l = tau
        tmax = np.max(tau_l)
        ind_m = np.where(tau_l == tmax)[0][0]
        tau_l = np.roll(tau_l, np.size(tau_l)/2- ind_m)
        if width != None:
            ldla = np.size(tau_l)/2 - width/2
            hdla = np.size(tau_l)/2 + width/2
        else:
            ldla = 0
            hdla = np.size(tau_l)
        plt.plot(np.arange(0,np.size(tau_l))*self.dvbin,np.exp(-tau_l))
        (low, high) = self._vel_width_bound(tau_l[ldla:hdla], tot_tau)
        low+=ldla
        high+=ldla
        if high - low > 0:
            plt.plot([low,low],[0,1])
            plt.plot([high,high],[0,1])
        if high - low > 150:
            tpos = low + 15
        else:
            tpos = high+15
        plt.text(tpos,0.5,r"$\delta v_{90} = "+str(np.round(high-low,1))+r"$")
        plt.ylim(-0.05,1.05)
        plt.xlim(low-70,high+70)
        return (low, high)

    def plot_spectrum_density_velocity(self, elem, ion, i):
        """Plot the spectrum of a line, centered on the deepest point,
           and marking the 90% velocity width."""
        #  Size of a single velocity bin
        tau = self.get_observer_tau(elem, ion,i)
        col_den = self.get_col_density(elem, ion)[i,:]
        #Deal with periodicity by making sure the deepest point is in the middle
        tau_l = np.ravel(tau)
        tmax = np.max(tau_l)
        ind_m = np.where(tau_l == tmax)[0][0]
        tau_l = np.roll(tau_l, np.size(tau_l)/2- ind_m)
        col_den = np.roll(col_den, np.size(tau_l)/2- ind_m)
        plt.subplot(211)
        (low, high) = self.plot_spectrum(tau)
        plt.xlim(low-50,high+50)
        plt.subplot(212)
        ind = np.where(col_den > 1)
        plt.semilogy(np.arange(0,np.size(tau_l))[ind]*self.dvbin,col_den[ind])
        plt.xlim(low-50,high+50)
        plt.yticks(np.array([1e10, 1e15,1e20]))

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

    def plot_sep_frac(self,elem = "Si", ion = 2, thresh = 1e-2, mindist = 15, dv = 0.1):
        """
        Plots the fraction of spectra in each velocity width bin which are separated.
        Threshold is as a percentage of the maximum value.
        mindist is in km/s
        """
        sep = self.get_separated(elem, ion, thresh,mindist)
        vels = self.vel_width(elem, ion)
        v_table = 10**np.arange(0, np.log10(np.max(vels)), dv)
        vbin = np.array([(v_table[i]+v_table[i+1])/2. for i in range(0,np.size(v_table)-1)])
        hist1 = np.histogram(vels, v_table)
        hist2 = np.histogram(vels[sep],v_table)
        hist1[0][np.where(hist1[0] == 0)] = 1
        plt.semilogx(vbin, hist2[0]/(1.*hist1[0]))

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
        ind = self.get_filt("Z", -1, None)
        self._plot_metallicity(met[ind],nbins,color,ls)

    def plot_species_metallicity(self, species, ion, nbins=20,color="blue", ls="-"):
        """Plot the distribution of metallicities from an ionic species"""
        met = self.get_ion_metallicity(species,ion)
        ind = self.get_filt(species, ion, None)
        self._plot_metallicity(met[ind],nbins,color,ls)

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

    def plot_Z_vs_vel_width(self,elem="Si", ion=2, color="blue"):
        """Plot the correlation between metallicity and velocity width"""
        ind = self.get_filt(elem, ion)
        vel = self.vel_width(elem, ion)[ind]
        met = self.get_metallicity()
        met = met[ind]
        #Ignore objects too faint to be seen or unresolved
        ind2 = np.where(np.logical_and(vel > 15, met > 1e-3))
        plt.loglog(vel[ind2],met[ind2], 'x',color=color)
        met = np.log10(met[ind2])
        vel = np.log10(vel[ind2])
        (intercept, slope, var) = lsq.leastsq(met,vel)
        print "sim corr: ",intercept, slope, np.sqrt(var)
        print "sim correlation: ",lsq.pearson(met, vel,intercept, slope)
        print "sim kstest: ",lsq.kstest(met, vel,intercept, slope)
        xx = np.logspace(np.min(met), np.max(met),15)
        plt.loglog(10**intercept*xx**slope, xx, color=color)
        plt.xlim(10,2e3)

    def plot_Z_vs_mass(self,color="blue"):
        """Plot the correlation between mass and metallicity, with a fit"""
        (halo, _) = self.find_nearest_halo()
        ind = np.where(halo > 0)
        halo = halo[ind]
        mass = self.sub_mass[halo]
        met = self.get_metallicity()[ind]
        plt.loglog(mass,met, 'x',color=color)
        met = np.log10(met)
        mass = np.log10(mass)
        mind = np.where(met > -4)
        met = met[mind]
        mass = mass[mind]
        (intercept, slope, var) = lsq.leastsq(met,mass)
        print "Z mass corr: ",intercept, slope, np.sqrt(var)
        print "Z mass correlation: ",lsq.pearson(met, mass,intercept, slope)
        print "Z mass kstest: ",lsq.kstest(met, mass,intercept, slope)
        xx = np.logspace(np.min(met), np.max(met),15)
        plt.loglog(10**intercept*xx**slope, xx, color="black", label=self.label)
        plt.ylim(1e-4,10)

    def plot_vel_vs_mass(self,elem, ion, color="blue"):
        """Plot the correlation between mass and metallicity, with a fit"""
        (halo, _) = self.find_nearest_halo()
        ind = np.where(halo > 0)
        halo = halo[ind]
        vel = self.vel_width(elem, ion)[ind]
        mass = self.sub_mass[halo]
        plt.loglog(mass,vel, 'x',color=color)
        vel = np.log10(vel)
        mass = np.log10(mass)
        (intercept, slope, var) = lsq.leastsq(vel,mass)
        print "Z vel corr: ",intercept, slope, np.sqrt(var)
        print "Z vel correlation: ",lsq.pearson(vel, mass,intercept, slope)
        print "Z vel kstest: ",lsq.kstest(vel, mass,intercept, slope)
        xx = np.logspace(np.min(vel), np.max(vel),15)
        plt.loglog(10**intercept*xx**slope, xx, color="black",label=self.label)

    def kstest(self, Zdata, veldata, elem="Si", ion=2):
        """Find the 2D KS test value of the Vel width and log metallicity with respect to an external dataset, veldata and Z data"""
        met = self.get_metallicity()
        ind = self.get_filt(elem, ion)
        met = np.log10(met[ind])
        vel = np.log10(self.vel_width(elem, ion)[ind])
        data2 = np.array([met,vel]).T
        data = np.array([np.log10(Zdata), np.log10(veldata)]).T
        return ks.ks_2d_2samp(data,data2)

    def plot_radius_vs_vel_width(self, elem="Si", ion=2, color="blue"):
        """Plot the velocity width vs the virial velocity of the hosting halo"""
        ind = self.get_filt(elem,ion)
        vel = self.vel_width(elem, ion)[ind]
        (halos, dists) = self.find_nearest_halo()
        radius = self.sub_radii[halos][ind]
        mass = self.sub_mass[halos][ind]
        ind2 = np.where(vel > 15)
        vel = vel[ind2]
        fromc = (dists[ind][ind2]/2.0/radius[ind2])**2
        fromc[np.where(fromc >= 0.95)] = 0.95
        virial = np.sqrt(4.302e-3*mass[ind2]/radius[ind2]/1000)
        plt.plot(virial, vel,'x',color=color)
        (intercept, slope, var) = lsq.leastsq(np.log10(virial),np.log10(vel))
        print "sim corr: ",intercept, slope, np.sqrt(var)
        print "sim correlation: ",lsq.pearson(np.log10(virial), np.log10(vel),intercept, slope)
        print "sim kstest: ",lsq.kstest(np.log10(virial), np.log10(vel),intercept, slope)
        xx = np.logspace(np.min(np.log10(virial)), np.max(np.log10(virial)),15)
        plt.loglog( xx,10**intercept*xx**slope, color=color)

        #(intercept, slope, var) = lsq.leastsq(virial,vel)
        #print "sim corr: ",intercept, slope, np.sqrt(var)
        #print "sim correlation: ",lsq.pearson(virial, vel,intercept, slope)
        #xx = np.linspace(np.min(virial), np.max(virial),15)
        #plt.loglog( xx,intercept+xx*slope, color=color)

    def plot_virial_vel_vs_vel_width(self,elem="Si", ion=2):
        """Plot a histogram of the velocity widths vs the halo virial velocity"""
        ind = self.get_filt(elem,ion)
        vel = self.vel_width(elem, ion)[ind]
        ind2 = np.where(vel > 15)
        vel = vel[ind2]
        (halos, _) = self.find_nearest_halo()
        #Grav constant 4.302e-3 parsec / solar mass (km/s)^2
        virial = np.sqrt(4.302e-3*self.sub_mass[halos][ind][ind2]/self.sub_radii[halos][ind][ind2]/1000)
        ind2 = np.where(vel < 300)
        (H, xedges) = np.histogram(np.log10(vel[ind2]/virial[ind2]), bins=20,normed=True)
        print "median v/vir: ",np.median(vel/virial)
        plt.semilogx(10**xedges[:-1], H, color="red")
        ind2 = np.where(vel > 300)
        (H, xedges) = np.histogram(np.log10(vel[ind2]/virial[ind2]), bins=20,normed=True)
        plt.semilogx(10**xedges[:-1], H, color="blue")


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

