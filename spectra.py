# -*- coding: utf-8 -*-
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
from _spectra_priv import _SPH_Interpolate, _near_lines,_Compute_Absorption

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
        self.KPC = 3.08568025e19
        MPC = self.KPC * 1000
        self.SIGMA_T = 6.652458558e-29
        self.PROTONMASS = 1.66053886e-27 # 1 a.m.u in kg
        self.SOLAR_MASS = 1.98892e30
        self.GAMMA = 5.0/3.0
        #Spectral data
        self.num = num
        self.base = base
        self.cofm = cofm
        self.axis = np.array(axis, dtype = np.int32)
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
        #Calculate omega_baryon (approximately)
        mass_dm = ff["Header"].attrs["MassTable"][1]*ff["Header"].attrs["NumPart_ThisFile"][1]
        mass_bar = np.sum(ff["PartType0"]["Masses"])
        self.omegab = mass_bar/(mass_bar+mass_dm)*self.OmegaM
        ff.close()

        # Conversion factors from internal units
        rscale = (self.KPC*self.atime)/self.hubble    # convert length to m
        mscale = (1.0e10*self.SOLAR_MASS)/self.hubble   # convert mass to kg
        self.dscale = mscale / rscale **3 # Convert density to kg / m^3
        #  Calculate the length scales to be used in the box
        self.Hz = 100.0*self.hubble * np.sqrt(self.OmegaM/self.atime**3 + self.OmegaLambda)
        self.vmax = self.box * self.Hz * rscale/ MPC # box size (kms^-1)
        self.dzgrid   = self.box * rscale / (1.*self.nbins) # bin size m
        self.dvbin = self.dzgrid * self.Hz / MPC # velocity bin size (kms^-1)
        #Species we can use
        self.species = ['H', 'He', 'C', 'N', 'O', 'Ne', 'Mg', 'Si', 'Fe']
        #Generate cloudy tables
        self.cloudy_table = convert_cloudy.CloudyTable(cloudy_dir, self.red)
        #Line data
        self.lines = line_data.LineData()
        #Empty dictionary to add results to
        self.metals = {}

    def SPH_Interpolate_metals(self, elem, ion, get_rho_H=False):
        """Interpolate particles to lines of sight, calculating density, temperature and velocity
        of various metal species along the line of sight.
        HI is special-cased.
        Note: the ionisation fraction is just cloudy. Some self-shielding might be useful.
        This is a wrapper which calls the C function.
        Arguments:
            elem - Element(s) to compute spectra of
            ion - Ion density to compute. Only one ion allowed right now
            get_rho_H - If this is true, compute the bare hydrogen density

        Returns:
            rho_H - hydrogen density along the line of sight if get_rho_H = True
            dictionary with a list of [density, velocity, temperature] for each species along the line of sight.
            Units are physical kg/m^3, km/s and K.
        """
        #Get array sizes
        (rho_H, rho_metal, vel_metal, temp_metal) =  self._interpolate_single_file(self.files[0], elem, ion, get_rho_H)
        #Do remaining files
        for fn in self.files[1:]:
            (trho_H, trho_metal, tvel_metal, ttemp_metal) =  self._interpolate_single_file(fn, elem, ion, get_rho_H)
            #Add new file
            if get_rho_H:
                rho_H += trho_H
            rho_metal += trho_metal
            vel_metal += tvel_metal
            temp_metal += ttemp_metal
            del trho_H
            del trho_metal
            del tvel_metal
            del ttemp_metal
        #Rescale units
        metals = self.rescale_units(rho_metal, vel_metal, temp_metal)
        if get_rho_H:
            rho_H *= self.dscale
            return [rho_H,]+ metals
        else:
            return metals


    def _interpolate_single_file(self,fn, elem, ion, rho_H):
        """Read arrays and perform interpolation for a single file"""
        nelem = self.species.index(elem)
        ff = h5py.File(fn)
        data = ff["PartType0"]
        pos = np.array(data["Coordinates"],dtype=np.float32)
        vel = np.array(data["Velocities"],dtype=np.float32)
        mass = np.array(data["Masses"],dtype=np.float32)
        u = np.array(data["InternalEnergy"],dtype=np.float32)
        ne = np.array(data["ElectronAbundance"],dtype=np.float32)
        hh = np.array(hsml.get_smooth_length(data),dtype=np.float32)
        #In kg/m^3
        den = np.array(data["Density"], dtype = np.float32)*self.dscale
        #In (hydrogen) atoms / cm^3
        den /= (self.PROTONMASS*100**3)
        #Find particles we care about
        ind = self.particles_near_lines(pos, hh)
        pos = pos[ind,:]
        vel = vel[ind,:]
        mass = mass[ind]
        u = u[ind]
        ne = ne[ind]
        hh = hh[ind]
        den = den[ind]
        #Get metallicity of this metal species
        try:
            metal_in = np.array(data["GFM_Metals"],dtype=np.float32)[:,nelem]
            #Deal with floating point roundoff - metal_in will sometimes be negative
            #10^-30 is Cloudy's definition of zero.
            metal_in[np.where(metal_in < 1e-30)] = 1e-30
            metal_in = metal_in[ind]
        except KeyError:
            #Some default abundances. H and He are primordial, the rest are Milky Way as given by wikipedia
            metal_abund = np.array([0.76, 0.24, 4.6e-3, 9.6e-4, 1.04e-2, 1.34e-3, 5.8e-4, 6.5e-4, 1.09e-3])
            metal_in = metal_abund[nelem]
        #Get density of this ion - we need to weight T and v by ion abundance
        #Cloudy density in physical H atoms / cm^3
        #Special case H1:
        if elem == 'H':
            if ion != 1:
                raise ValueError
            # Hydrogen mass frac in the data array
            metal_in = np.array(data["NeutralHydrogenAbundance"],dtype=np.float32)[ind]*metal_in
        else:
            metal_in = self.cloudy_table.ion(elem, ion, metal_in, den)*metal_in
        ff.close()
        if rho_H:
            return _SPH_Interpolate(1,self.nbins, self.box, pos, vel, mass, u, ne, metal_in, hh, self.axis, self.cofm)
        else:
            return (None,)+_SPH_Interpolate(0,self.nbins, self.box, pos, vel, mass, u, ne, metal_in, hh, self.axis, self.cofm)

    def particles_near_lines(self, pos, hh):
        """Filter a particle list, returning an index list of those near sightlines"""
        ind = _near_lines(self.box, pos, hh, self.axis, self.cofm)
        return ind

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

    def compute_absorption(self,elem, ion, ll, rho, vel, temp):
        """Computes the absorption spectrum (tau (u) ) from a binned set of interpolated
        densities, velocities and temperatures.

        Optical depth is given by:
        tau (u) = sigma_X c / H(z) int_infty^infty n_x(x) V( u - x - v_pec, b(x) ) dx
        where V is the Voigt profile, b(x)^2 = 2k_B T /m_x c^2 is the velocity dispersion.
        and v_pec is the peculiar velocity.
        sigma_X is the cross-section for this transition.
        """
        #Get line data
        line = self.lines[(elem,ion)][ll]
        mass = self.lines.get_mass(elem)
        #Don't forget to convert line width from A to m!
        tau = _Compute_Absorption(rho, vel, temp, self.nbins, self.Hz, self.hubble, self.box, self.atime,line.lambda_X*1e-10, line.gamma_X, line.fosc_X,mass)
        return tau

    def find_max_tau(self, elem, ion, rho, vel, temp):
        """Find which of the transitions gives the largest maximal optical depth."""
        line = self.lines[(elem,ion)]
        mass = self.lines.get_mass(elem)
        maxes = [np.max(_Compute_Absorption(rho, vel, temp, self.nbins, self.Hz, self.hubble, self.box, self.atime,ll.lambda_X*1e-10, ll.gamma_X, ll.fosc_X,mass)) for ll in line]
        return np.where(maxes == np.max(maxes))

    def compute_absorption_python(self,elem, ion, ll, rho, vel, temp):
        """Computes the absorption spectrum (tau (u) ) from a binned set of interpolated
        densities, velocities and temperatures.

        Optical depth is given by:
        tau (u) = sigma_X c / H(z) int_infty^infty n_x(x) V( u - x - v_pec, b(x) ) dx
        where V is the Voigt profile, b(x)^2 = 2k_B T /m_x c^2 is the velocity dispersion.
        and v_pec is the peculiar velocity.
        sigma_X is the cross-section for this transition.
        """
        #Get line data
        line = self.lines[(elem,ion)][ll]
        line.lambda_X*=1e-10
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

    def get_tau(self, elem, ion, ll):
        """Get the optical depth for a particular element out of:
           (He, C, N, O, Ne, Mg, Si, Fe)
           and some ion number
        """
        #generate metal and hydrogen spectral densities
        #Indexing is: rho_metals [ NSPECTRA, NBIN ]
        [rho, vel, temp] = self.SPH_Interpolate_metals(elem, ion)

        #Compute tau for this metal ion
        (nlos, nbins) = np.shape(rho)
        tau = np.array([self.compute_absorption(elem, ion, ll, rho[n,:], vel[n,:], temp[n,:]) for n in xrange(0, nlos)])
        self.metals[(elem, ion)] = [rho, vel, temp, tau]
        return tau

    def vel_width(self, tau):
        """Find the velocity width of a line
           defined as the width of 90% of the integrated optical depth.
           This is a little complicated by periodic boxes,
           so we internally cycle the line until the deepest absorption
           is in the middle"""
        #  Size of a single velocity bin
        tot_tau = np.sum(tau,axis = 1)
        vel_width = np.zeros(np.shape(tot_tau))
        for ll in np.arange(0, np.shape(tau)[0]):
            #Deal with periodicity by making sure the deepest point is in the middle
            tau_l = tau[ll,:]
            max = np.max(tau_l)
            ind_m = np.where(tau_l == max)[0][0]
            tau_l = np.roll(tau_l, np.size(tau_l)/2- ind_m)
            cum_tau = np.cumsum(tau_l)
            ind_low = np.where(cum_tau > 0.05 * tot_tau[ll])
            low = ind_low[0][0]
            ind_high = np.where(cum_tau > 0.95 * tot_tau[ll])
            high = ind_high[0][0]
            vel_width[ll] = self.dvbin*(high-low)
        #Return the width
        return vel_width

    def NHI(self):
        """Get the neutral hydrogen column densities for each line"""
        [rho, vel, temp] = self.SPH_Interpolate_metals('H', 1)
        #Column density in kg / m^3 (* box [comoving kpc/h in physical m])
        col_rho = np.sum(rho, axis=1)*self.box/(1.*self.nbins)*self.KPC*(1+self.red)/self.hubble
        #In atoms / m^2
        col_rho /= self.PROTONMASS
        #In atoms / cm^2
        col_rho /= 100**2
        return col_rho

