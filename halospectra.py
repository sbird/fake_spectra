# vim: set fileencoding=utf-8
"""Class to gather and analyse various metal line statistics"""

import numpy as np
import hdfsim
import halocat
import spectra
import matplotlib.pyplot as plt

class HaloSpectra(spectra.Spectra):
    """Generate metal line spectra from simulation snapshot"""
    def __init__(self,num, base, minpart = 400, nbins = 1024, cloudy_dir="/home/spb/codes/ArepoCoolingTables/tmp_spb/"):
        #Load halo centers to push lines through them
        f = hdfsim.get_file(num, base, 0)
        self.OmegaM = f["Header"].attrs["Omega0"]
        self.box = f["Header"].attrs["BoxSize"]
        self.npart=f["Header"].attrs["NumPart_Total"]+2**32*f["Header"].attrs["NumPart_Total_HighWord"]
        min_mass = self.min_halo_mass(minpart)
        f.close()
        (ind, self.sub_mass, cofm, self.sub_radii) = halocat.find_wanted_halos(num, base, min_mass)
        self.NumLos = np.size(self.sub_mass)
        #Random integers from [1,2,3]
#         axis = np.random.random_integers(3, size = self.NumLos)
        #All through y axis
        axis = 2*np.ones(self.NumLos)
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

    def vel_width_hist(self,tau, dv=0.2):
        """
        This computes the DLA column density function, which is the number
        of absorbers per sight line with velocities in the interval
        [v, v+dv] at the absorption distance X.
        Absorption distance is simply a single simulation box.
        A sightline is assumed to be equivalent to one grid cell.
        That is, there is presumed to be only one halo in along the sightline
        encountering a given halo.

        So we have f(N) = d n/ dv dX
        and n(N) = number of absorbers per sightline in this velocity bin.
        ie, f(N) = n / Δv / ΔX
        Note f(N) has dimensions of km/s, because v has units of km/s and X is dimensionless.

        Parameters:
            tau - optical depth along sightline
            dv - bin spacing

        Returns:
            (v, f_table) - v (binned in log) and corresponding f(N)
        """
        vel_width = self.vel_width(tau)
        nlos = np.shape(tau)[0]
        v_table = np.arange(0, np.log10(np.max(vel_width)), dv)
        bin = np.array([(v_table[i]+v_table[i+1])/2. for i in range(0,np.size(v_table)-1)])
        dX=self.absorption_distance()
        nn = np.histogram(np.log10(vel_width),v_table)[0] / (1.*nlos)
        vels=nn/(dv*10**bin*dX)
        return (10**bin, vels)

    def plot_vel_width(self, tau):
        (bin, vels) = self.vel_width_hist(tau)
        plt.loglog(bin, vels)
