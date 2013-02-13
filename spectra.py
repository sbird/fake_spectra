"""Python module for generating fake spectra from an N-body catalogue.

Note in Arepo we have GFM_Metals and GFM_Metallicity.

GFM_Metallicity is the total mass in species not H or He
per unit gas mass (and is used for cooling).

GFM_Metals is a 9-component array of species:
H, He, C, N, O, Ne, Mg, Si, Fe

Because these are not all the species, GFM_Metals will not sum to 1
and sum(GFM_Metals[2:])  < GFM_Metallicity

However, it should be true that
1- sum(GFM_Metals[:]) + GFM_Metals[2:] = GFM_Metallicity

Also note that there is some instability at very low metallicities - the code will often return +-1e-20.
"""


import numpy as np
import hsml
import math
import convert_cloudy
import line_data
import h5py
import hdfsim
from _spectra_priv import _SPH_Interpolate

class Spectra:
    """Class to interpolate particle densities along a line of sight and calculate their absorption
        Arguments:
            num - Snapshot number
            base - Name of base directory for snapshot
	        cofm - table of los positions, as [n, 3] shape array.
            axis - axis along which to put the sightline
	        nbins (optional) - number of bins in each spectrum
            cloudy_dir (optional) - Directory containing cloudy tables for ionisation fraction calculations"""
    def __init__(self,num, base,cofm, axis, nbins=1024, cloudy_dir="/home/spb/codes/ArepoCoolingTables/tmp_spb/"):
        #Various physical constants
        #Speed of light
        self.light = 2.99e8
        #Boltzmann constant
        self.BOLTZMANN = 1.3806504e-23
        KPC = 3.08568025e19
        MPC = KPC * 1000
        self.SIGMA_T = 6.652458558e-29
        self.PROTONMASS = 1.66053886e-27 # 1 a.m.u in kg
        self.SOLAR_MASS = 1.98892e30
        self.GAMMA = 5.0/3.0
        #Spectral data
        self.num = num
        self.base = base
        self.cofm = cofm
        self.axis = axis
        self.nbins = nbins
        #Snapshot data
        self.files = hdfsim.get_all_files(num, base)
        ff = h5py.File(self.files[0])
        self.box = ff["Header"].attrs["BoxSize"]
        self.red = ff["Header"].attrs["Redshift"]
        self.atime = ff["Header"].attrs["Time"]
        self.hubble = ff["Header"].attrs["HubbleParam"]
        self.OmegaM = ff["Header"].attrs["Omega0"]
        self.OmegaLambda = ff["Header"].attrs["OmegaLambda"]
        ff.close()
        # Conversion factors from internal units
        rscale = (KPC*self.atime)/self.hubble    # convert length to m
        mscale = (1.0e10*self.SOLAR_MASS)/self.hubble   # convert mass to kg
        self.dscale = mscale / rscale **3 # Convert density to kg / m^3
        #  Calculate the length scales to be used in the box
        Hz = 100.0*self.hubble * np.sqrt(self.OmegaM/self.atime**3 + self.OmegaLambda)
        self.vmax = self.box * Hz * rscale/ MPC # box size (kms^-1)
        self.dzgrid   = self.box * rscale / (1.*self.nbins) # bin size m
        self.dvbin = self.dzgrid * Hz / MPC # velocity bin size (kms^-1)
        #Species we can use
        self.species = ['H', 'He', 'C', 'N', 'O', 'Ne', 'Mg', 'Si', 'Fe']
        #Generate cloudy tables
        self.cloudy_table = convert_cloudy.CloudyTable(cloudy_dir)
        #Line data
        self.lines = line_data.LineData()

    def SPH_Interpolate_metals(self, elem, ion):
        """Interpolate particles to lines of sight, calculating density, temperature and velocity
        of various metal species along the line of sight.

        This is a wrapper which calls the C function.
        Arguments:
            num - Snapshot number
            base - Name of base directory for snapshot
	        cofm - table of los positions, as [n, 3] shape array.
            axis - axis along which to put the sightline
	        nbins - number of bins in each spectrum
            elem - Element to compute spectra of
            ion - Ion density to compute
            cloudy_table - Object containing cloudy tables for ionisation fraction calculations

        Returns:
            rho_H - hydrogen density along the line of sight
            dictionary with a list of [density, velocity, temperature] for each species along the line of sight.
        """
        #Get array sizes
        (rho_H, rho_metal, vel_metal, temp_metal) =  self._interpolate_single_file(self.files[0], elem, ion)
        #Do remaining files
        for fn in self.files[1:]:
            (trho_H, trho_metal, tvel_metal, ttemp_metal) =  self._interpolate_single_file(fn, elem, ion)
            #Add new file
            rho_H += trho_H
            rho_metal += trho_metal
            vel_metal += tvel_metal
            temp_metal += ttemp_metal
            del trho_H
            del trho_metal
            del tvel_metal
            del ttemp_metal
        metals = {}
        #Rescale units
        rho_H *= self.dscale
        for mm in np.arange(0,np.shape(rho_metal)[2]):
            metals[elem] = [self.rescale_units(rho_metal[:,:,mm], vel_metal[:,:,mm], temp_metal[:,:,mm])]
        return (rho_H, metals)


    def _interpolate_single_file(self,fn, elem, ion):
        """Read arrays and perform interpolation for a single file"""
        nelem = [self.species.index(ee) for ee in elem]
        ff = h5py.File(fn)
        data = ff["PartType0"]
        pos = np.array(data["Coordinates"],dtype=np.float32)
        vel = np.array(data["Velocities"],dtype=np.float32)
        mass = np.array(data["Masses"],dtype=np.float32)
        u = np.array(data["InternalEnergy"],dtype=np.float32)
        ne = np.array(data["ElectronAbundance"],dtype=np.float32)
        hh = np.array(hsml.get_smooth_length(data),dtype=np.float32)
        #Get metallicity of this metal species
        metal_in = np.array(data["GFM_Metals"],dtype=np.float32)[:,nelem]
        #In kg/m^3
        den = np.array(data["Density"], dtype = np.float32)*self.dscale
        #In (hydrogen) atoms / cm^3
        den /= (self.PROTONMASS*100**3)
        ff.close()
        #Deal with floating point roundoff - metal_in will sometimes be negative
        metal_in[np.where(metal_in < 0)] = 0
        #Get density of this ion - we need to weight T and v by ion abundance
        #Cloudy density in physical H atoms / cm^3
        ion = self.cloudy_table.ion(elem, ion, self.red, metal_in, den)
        return _SPH_Interpolate(self.nbins, self.box, pos, vel, mass, u, ne, metal_in*ion, hh, self.axis, self.cofm)

    def rescale_units(self, rho, vel, temp):
        """Rescale the units of the arrays from internal gadget units to
        physical kg/m^3, km/s and K.
            Only do this ONCE."""
        # Conversion factors from internal units
        vscale = np.sqrt(self.atime)        #convert velocity to kms^-1
        # Rescale density and vel. temp is already in K
        # vel and temp are calculated weighted by density. Undo this.
        ind = np.where(rho > 0)
        vel[ind] *= vscale/rho[ind]
        temp[ind] /= rho[ind]
        rho[ind] *= self.dscale
        #If there are no particles in this bin, rho will be zero.
        #In this case, we set temp and veloc arbitrarily to one,
        #to avoid nans propagating. Zero rho will imply zero absorption
        #anyway.
        ind = np.where(rho == 0)
        vel[ind]=1
        temp[ind]=1
        return [rho, vel, temp]

    def compute_absorption(self,elem, ion, rho, vel, temp):
        """Computes the absorption spectrum (tau (u) ) from a binned set of interpolated
        densities, velocities and temperatures.

        Optical depth is given by:
        tau (u) = sigma_X c / H(z) int_infty^infty n_x(x) V( u - x - v_pec, b(x) ) dx
        where V is the Voigt profile, b(x)^2 = 2k_B T /m_x c^2 is the velocity dispersion.
        and v_pec is the peculiar velocity.
        sigma_X is the cross-section for this transition.
        """
        #Get line data
        line = self.lines.get_line(elem,ion)
        mass = self.lines.get_mass(elem)
        tau = np.zeros(self.nbins)

        #Absorption cross-sections m^2
        sigma_X  = np.sqrt(3.0*math.pi*self.SIGMA_T/8.0) * line.lambda_X  * line.fosc_X
        # Prefactor for optical depth
        A_H1 = sigma_X*self.light*self.dzgrid/np.sqrt(math.pi)
        #Compute the spectra optical depth
        for i in xrange(0, self.nbins):
            uu = self.dvbin*1.e3*np.arange(0,self.nbins)
            uu += vel*1.e3
            # Note this is indexed with i, above with j!
            # This is the difference in velocities between two clouds
            # on the same sightline
            vdiff  = np.abs(self.dvbin*i*1.0e3 - uu)  # ms^-1
            ind = np.where(vdiff > self.vmax *1.e3 /2.)
            vdiff[ind] = self.vmax*1.e3 - vdiff[ind]
            #Impact parameter
            bb = np.sqrt(2.0*self.BOLTZMANN*temp/(mass*self.PROTONMASS))
            T0 = (vdiff/bb)**2
            T1 = np.exp(-T0)
            aa_H1 = line.gamma_X*line.lambda_X/(4.0*math.pi*bb)
            T2 = 1.5/T0
            ind = np.where(T0 > 1.e-6)
            profile = np.array(T1)
            # Voigt profile: Tepper-Garcia, 2006, MNRAS, 369, 2025
            #This appears to break down for high T.
            profile[ind] = T1[ind] - aa_H1[ind]/np.sqrt(math.pi)/T0[ind]*(T1[ind]**2*(4.0*T0[ind]**2 + 7.0*T0[ind] + 4.0 + T2[ind]) - T2[ind] -1.0)
            tau[i] = np.sum(A_H1  * rho  * profile /(mass*self.PROTONMASS*bb))

        return tau

    def get_tau(self, elem, ion):
        """Get the optical depth for a particular element out of:
           (He, C, N, O, Ne, Mg, Si, Fe)
           and some ion number
           NOTE: May wish to special-case SiIII at some point
        """
        #generate metal and hydrogen spectral densities
        #Indexing is: rho_metals [ NSPECTRA, NBIN ]
        (rho_H, metals) = self.SPH_Interpolate_metals(elem, ion)

        species = metals[elem]
        #Compute tau for this metal ion
        (nlos, nbins) = np.shape(species[0])
        tau_metal=np.empty((nlos, nbins))
        for n in xrange(0,nlos):
            tau_metal[n,:] = self.compute_absorption(species[0][n,:], species[1][n,:], species[2][n,:],elem, ion)
        return tau_metal

    def vel_width(self, tau):
        """Find the velocity width of a line"""
        #  Size of a single velocity bin
        tot_tau = np.sum(tau,axis = 1)
        cum_tau = np.cumsum(tau,axis = 1)
        vel_width = np.zeros(np.shape(tot_tau))
        for ll in np.arange(0, np.shape(tau)[1]):
            ind_low = np.where(cum_tau[ll,:] > 0.05 * tot_tau[ll])
            low = ind_low[0][-1]
            ind_high = np.where(cum_tau[ll,:] > 0.95 * tot_tau[ll])
            high = ind_high[0][0]
            vel_width[ll] = self.dvbin*(high-low)
        #Return the width
        return vel_width


