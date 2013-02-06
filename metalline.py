"""Class to gather and analse various metal line statistics"""

import numpy as np
import hdfsim
import convert_cloudy
import line_data
import spectra

class MetalLines:
    """Generate metal line spectra from simulation snapshot"""
    def __init__(self,num, base, los_table, cloudy_dir="/home/spb/codes/ArepoCoolingTables/tmp_spb/", nbins = 1024):
        f = hdfsim.get_file(num, base, 0)
        self.box = f["Header"].attrs["BoxSize"]
        self.hubble = f["Header"].attrs["HubbleParam"]
        self.atime = f["Header"].attrs["Time"]
        self.redshift = f["Header"].attrs["Redshift"]
        Omega0 = f["Header"].attrs["Omega0"]
        OmegaLambda = f["Header"].attrs["OmegaLambda"]
        f.close()
        self.Hz = 100.0*self.hubble * np.sqrt(Omega0/self.atime**3 + OmegaLambda)
        self.nbins = nbins
        self.xbins = np.arange(0,self.nbins)*self.box/self.nbins
        self.species = ['He', 'C', 'N', 'O', 'Ne', 'Mg', 'Si', 'Fe']
        self.NumLos = np.size(los_table.axis)
        #Line data
        self.lines = line_data.LineData()
        #generate metal and hydrogen spectral densities
        #Indexing is: rho_metals [ NSPECTRA, NBIN ]
        (self.rho_H, self.metals) = spectra.SPH_Interpolate_metals(num, base, los_table, nbins)
        #rescale H density
        self.rho_H = spectra.rescale_units_rho_H(self.rho_H, self.hubble, self.atime)
        #Rescale metals
        for (key, value) in self.metals.iteritems():
            mass = self.lines.get_mass(key)
            value.rescale_units(self.hubble, self.atime, mass)

        #Generate cloudy tables
        self.cloudy = convert_cloudy.CloudyTable(cloudy_dir)


    def get_lines(self, elem, ion):
        """Get the optical depth for a particular element out of:
           (He, C, N, O, Ne, Mg, Si, Fe)
           and some ion number
           NOTE: May wish to special-case SiIII at some point
        """
        species = self.metals[elem]
        line = self.lines.get_line(elem,ion)
        mass = self.lines.get_mass(elem)
        #To convert H density from kg/m^3 to atoms/cm^3
        conv = 1./(spectra.PROTONMASS*100**3)
        ion_density = np.array(species.rho)
        #Compute tau for this metal ion
        tau_metal=np.empty(np.shape(species.rho))
        for n in np.arange(0,self.NumLos):
            #For the density parameter use the hydrogen density at this pixel
            #For metallicity pass the metallicity of this species at this bin (rho_Z/ rho_H) and it will be converted to cloudy format
            ind = np.where(ion_density[n,:] > 0)
            for i in np.array(ind).ravel():
                ion_density[n,i] *= self.cloudy.ion(elem, ion, self.redshift, ion_density[n,i]/self.rho_H[n,i], self.rho_H[n,i]/conv)
            tau_metal[n] = spectra.compute_absorption(self.xbins, ion_density[n], species.vel[n], species.temp[n],line,self.Hz,self.hubble, self.box, self.atime,mass)
        return tau_metal

    def vel_width(self, tau):
        """Find the velocity width of a line"""
        #  Size of a single velocity bin
        dvbin = self.box / (1.*self.nbins) * self.Hz *self.atime /self.hubble / 1000 # velocity bin size (kms^-1)

        tot_tau = np.sum(tau)
        cent = np.where(tau == np.max(tau))
        tot = tau[cent]
        if tau[cent] > 0.9 *tot_tau:
            return dvbin
        i = 0
        #Extend the region of interest until we have enough tau
        while tot < 0.9*tot_tau:
            i+=1
            tot += tau[cent-i] + tau[cent+i]
        #Return the width
        return dvbin*i*2

