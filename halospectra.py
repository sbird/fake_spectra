# -*- coding: utf-8 -*-
"""Class to gather and analyse various metal line statistics"""

import numpy as np
import hdfsim
import math
import halocat
import spectra
import matplotlib.pyplot as plt

class HaloSpectra(spectra.Spectra):
    """Generate metal line spectra from simulation snapshot"""
    def __init__(self,num, base, repeat = 3, minpart = 400, nbins = 1024, cloudy_dir="/home/spb/codes/ArepoCoolingTables/tmp_spb/"):
        #Load halos to push lines through them
        f = hdfsim.get_file(num, base, 0)
        self.OmegaM = f["Header"].attrs["Omega0"]
        self.box = f["Header"].attrs["BoxSize"]
        self.npart=f["Header"].attrs["NumPart_Total"]+2**32*f["Header"].attrs["NumPart_Total_HighWord"]
        min_mass = self.min_halo_mass(minpart)
        f.close()
        (ind, self.sub_mass, cofm, self.sub_radii) = halocat.find_wanted_halos(num, base, min_mass)
        self.NumLos = np.size(self.sub_mass)*repeat
        #All through y axis
        axis = np.ones(self.NumLos)
        axis[self.NumLos/3:2*self.NumLos/3] = 2
        axis[2*self.NumLos/3:self.NumLos] = 3
        cofm = np.repeat(cofm,repeat,axis=0)
        axis = np.repeat(axis,repeat/3)
        #Perturb the sightlines within a sphere of half the virial radius.
        #We want a representative sample of DLAs.
        maxr = self.sub_radii/2.
        #Generate random sphericals
        theta = 2*math.pi*np.random.random_sample(self.NumLos)-math.pi
        phi = 2*math.pi*np.random.random_sample(self.NumLos)
        rr = np.repeat(maxr,repeat)*np.random.random_sample(self.NumLos)
        #Add them to halo centers
        cofm[:,0]+=rr*np.sin(theta)*np.cos(phi)
        cofm[:,1]+=rr*np.sin(theta)*np.sin(phi)
        cofm[:,2]+=rr*np.cos(theta)
        #Re-seed for repeatability
        np.random.seed(23)

        spectra.Spectra.__init__(self,num, base, cofm, axis, nbins, cloudy_dir)

    def min_halo_mass(self, minpart = 400):
        """Min resolved halo mass in internal Gadget units (1e10 M_sun)"""
        #This is rho_c in units of h^-1 1e10 M_sun (kpc/h)^-3
        rhom = 2.78e+11* self.OmegaM / 1e10 / (1e3**3)
        #Mass of an SPH particle, in units of 1e10 M_sun, x omega_m/ omega_b.
        target_mass = self.box**3 * rhom / self.npart[0]
        min_mass = target_mass * minpart
        return min_mass

    def absorption_distance(self):
        """Compute X(z), the absorption distance per sightline (eq. 9 of Nagamine et al 2003)
        in dimensionless units."""
        #h * 100 km/s/Mpc in h/s
        h100=3.2407789e-18
        #Units: h/s   s/m                        kpc/h      m/kpc
        return h100/self.light*(1+self.red)**2*self.box*self.KPC

    def absorption_distance_dz(self):
        """Compute dX/dz = H_0 (1+z)^2 / H(z) (which is independent of h)"""
        zp1 = 1+self.red
        return zp1**2/np.sqrt(self.OmegaM*zp1**3+(1-self.OmegaM))

    def vel_width_hist(self, tau, dv=0.1, col_rho=None):
        """
        Compute a histogram of the velocity widths of our spectra, with the purpose of
        comparing to the data of Prochaska 2008.

        Note this does not match Pontzen 2008, who multiply by the DLA fraction (0.065) obtained from the cddf.

        So we have f(N) = d n/ dv
        and n(N) = number of absorbers per sightline in this velocity bin.
        Note f(N) has dimensions of s/km, because v has units of km/s.

        Parameters:
            tau - optical depth along sightline
            dv - bin spacing
            col_rho - Prochaska used a subsample of spectra containing a DLA.
                      If this value is not None, use it as HI column density measurements to select
                      such a subsample. This does not change the resulting velocity widths.

        Returns:
            (v, f_table) - v (binned in log) and corresponding f(N)
        """
        #Remember this is not in log...
        if col_rho != None:
          ind = np.where(col_rho > 10**20.3)
          tau = tau[ind]
        vel_width = self.vel_width(tau)
        nlos = np.shape(vel_width)[0]
        print 'nlos = ',nlos
        v_table = 10**np.arange(0, np.log10(np.max(vel_width)), dv)
        bin = np.array([(v_table[i]+v_table[i+1])/2. for i in range(0,np.size(v_table)-1)])
        nn = np.histogram(vel_width,v_table)[0] / (1.*nlos)
        width = np.array([v_table[i+1]-v_table[i] for i in range(0,np.size(v_table)-1)])
        vels=nn/width
        return (bin, vels)

    def plot_vel_width(self, tau, dv=0.1, col_rho=None):
        """Plot the velocity widths of this snapshot
           Parameters:
            tau - optical depth along sightline
            dv - bin spacing
            col_rho - Prochaska used a subsample of spectra containing a DLA.
                      If this value is not None, use it as HI column density measurements to select
                      such a subsample. This does not change the resulting velocity widths.

        """
        (bin, vels) = self.vel_width_hist(tau, dv, col_rho)
        plt.semilogy(bin, vels)

    def get_col_density(self, elem, ion):
        """Get the column density in each pixel for a given species, assuming rho was already calculated"""
        rho = self.metals[(elem, ion)][0]
        #Size of a bin in physical m
        binsz = self.box/(1.*self.nbins)*self.KPC*(1+self.red)/self.hubble
        #Convert from physical kg/m^2 to atoms/cm^2
        convert = 1./self.PROTONMASS/1e4/self.lines.get_mass(elem)
        return rho*binsz*convert

    def plot_spectrum(self, tau, i):
        """Plot the spectrum of a line, centered on the deepest point,
           and marking the 90% velocity width."""
        #  Size of a single velocity bin
        tot_tau = np.sum(tau[i,:])
        #Deal with periodicity by making sure the deepest point is in the middle
        tau_l = tau[i,:]
        max = np.max(tau_l)
        ind_m = np.where(tau_l == max)[0][0]
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
