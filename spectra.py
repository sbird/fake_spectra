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
import h5py
import convert_cloudy
import line_data
from _spectra_priv import _SPH_Interpolate

#Various physical constants
#Speed of light
C = 2.99e8
#Boltzmann constant
BOLTZMANN = 1.3806504e-23
KPC = 3.08568025e19
MPC = KPC * 1000
SIGMA_T = 6.652458558e-29
PROTONMASS = 1.66053886e-27 # 1 a.m.u
SOLAR_MASS = 1.98892e30
GAMMA = 5.0/3.0

def SPH_Interpolate_metals(data, los_table, nbins, box):
    """Interpolate particles to lines of sight, calculating density, temperature and velocity
    of various metal species along the line of sight.

    This is a wrapper which calls the C function.
    Arguments:
    	data - HDF5 dataset from snapshot. Use f["PartType0"]
	    los_table - table of los positions. should have member arrays x, y, z and axis.
	    nbins - number of bins in each spectrum
	    box - box size
        hub - hubble constant (eg 0.71)
        atime - 1/(1+z)

    Returns:
        rho_H - hydrogen density along the line of sight
        dictionary of Species classes, specifying density, temperature and velocity of the metal species along the line of sight.
    """
    pos = np.array(data["Coordinates"],dtype=np.float32)
    vel = np.array(data["Velocities"],dtype=np.float32)
    mass = np.array(data["Masses"],dtype=np.float32)
    u = np.array(data["InternalEnergy"],dtype=np.float32)
    ne = np.array(data["ElectronAbundance"],dtype=np.float32)
    hh = np.array(hsml.get_smooth_length(data),dtype=np.float32)
    xx=np.array(los_table.xx, dtype=np.float32)
    yy=np.array(los_table.yy, dtype=np.float32)
    zz=np.array(los_table.zz, dtype=np.float32)
    axis=np.array(los_table.axis, dtype=np.int32)
    #We exclude hydrogen
    metal_in = np.array(data["GFM_Metals"],dtype=np.float32)[:,1:]
    #Deal with floating point roundoff - metal_in will sometimes be negative
    metal_in[np.where(np.abs(metal_in) < 1e-10)] = 0
    (rho_H, rho_metal, vel_metal, temp_metal) =  _SPH_Interpolate(nbins, box, pos, vel, mass, u, ne, metal_in, hh, axis, xx, yy, zz)
    species = ['He', 'C', 'N', 'O', 'Ne', 'Mg', 'Si', 'Fe']
    metals = {}
    for mm in np.arange(1,np.shape(metal_in)[1]):
        metals[species[mm]] = Species(rho_metal[:,:,mm], vel_metal[:,:,mm], temp_metal[:,:,mm])
    return (rho_H, metals)

class Species:
    """Convenience class to aggregate rho, vel and temp for a class"""
    def __init__(self, rho, vel, temp):
        self.rho = rho
        self.vel = vel
        self.temp = temp

    def rescale_units(self,h100, atime, mass):
        """Rescale the units of the arrays from internal gadget units to
        physical kg/m^3, km/s and K.
        Arguments:
            h100 = hubble constant
            atime = scale factor
            mass = mass of this species in amu
            Needed for the conversion between comoving kpc/h to physical m.
            Only do this ONCE."""
        # Conversion factors from internal units
        rscale = (KPC*atime)/h100    # convert length to m
        vscale = np.sqrt(atime)        #convert velocity to kms^-1
        mscale = (1.0e10*SOLAR_MASS)/h100   # convert mass to kg
        escale = 1.0e6           # convert energy/unit mass to J kg^-1
        # convert (with mu) T to K
        tscale = ((GAMMA-1.0) * mass * PROTONMASS * escale ) / BOLTZMANN
        # Rescale density, vel and temp
        # vel and temp are calculated weighted by density. Undo this.
        ind = np.where(self.rho > 0)
        self.vel[ind] *= vscale/self.rho[ind]
        self.temp[ind] *= tscale/self.rho[ind]
        self.rho[ind] *= mscale/rscale**3
        #If there are no particles in this bin, rho will be zero.
        #In this case, we set temp and veloc arbitrarily to one,
        #to avoid nans propagating. Zero rho will imply zero absorption
        #anyway.
        ind = np.where(self.rho == 0)
        self.vel[ind]=1
        self.temp[ind]=1

def rescale_units_rho_H(rho_H, h100, atime):
    """Rescale the units on hydrogen density"""
    # Conversion factors from internal units
    rscale = (KPC*atime)/h100    # convert length to m
    mscale = (1.0e10*SOLAR_MASS)/h100   # convert mass to kg
    return rho_H * mscale/rscale**3

def compute_absorption(xbins, rho, vel, temp, line, Hz, h100, box100, atime, mass):
    """Computes the absorption spectrum (tau (u) ) from a binned set of interpolated
    densities, velocities and temperatures.
    xbins are the positions of each bin along the sightline. A good default is
    xbins = np.range(0,NBINS)*Box/NBINS

    Optical depth is given by:
    tau (u) = sigma_X c / H(z) int_infty^infty n_x(x) V( u - x - v_pec, b(x) ) dx
    where V is the Voigt profile, b(x)^2 = 2k_B T /m_x c^2 is the velocity dispersion.
    and v_pec is the peculiar velocity.
    sigma_X is the cross-section for this transition.
    """
    nbins = np.size(xbins)
    tau = np.zeros(nbins)
    #  Conversion factors from internal units
    rscale = (KPC*atime)/h100  # convert length to m
    #  Calculate the length scales to be used in the box
    vmax = box100 * Hz * rscale/ MPC # box size (kms^-1)
    dzgrid   = box100 * rscale / 1.*nbins # bin size m
    dvbin = dzgrid * Hz / MPC # velocity bin size (kms^-1)

    #Absorption cross-sections m^2
    sigma_X  = np.sqrt(3.0*math.pi*SIGMA_T/8.0) * line.lambda_X  * line.fosc_X
    # Prefactor for optical depth
    A_H1 = sigma_X*C*dzgrid/np.sqrt(math.pi)
    #Compute the spectra optical depth
    for i in np.arange(0, nbins):
        uu = dvbin*1.e3*np.arange(0,nbins)
        uu += vel*1.e3
        # Note this is indexed with i, above with j!
        # This is the difference in velocities between two clouds
        # on the same sightline
        vdiff  = np.abs(dvbin*i*1.0e3 - uu)  # ms^-1
        ind = np.where(vdiff > vmax *1.e3 /2.)
        vdiff[ind] = vmax*1.e3 - vdiff[ind]
        #Impact parameter
        bb = np.sqrt(2.0*BOLTZMANN*temp/(mass*PROTONMASS))
        T0 = (vdiff/bb)**2
        T1 = np.exp(-T0)
        # Voigt profile: Tepper-Garcia, 2006, MNRAS, 369, 2025
        aa_H1 = line.gamma_X*line.lambda_X/(4.0*math.pi*bb)
        T2 = 1.5/T0
        ind = np.where(T0 > 1.e-6)
        profile = np.array(T1)
        profile[ind] = T1[ind] - aa_H1[ind]/np.sqrt(math.pi)/T0[ind]*(T1[ind]**2*(4.0*T0[ind]**2 + 7.0*T0[ind] + 4.0 + T2[ind]) - T2[ind] -1.0)
        tau[i] += np.sum(A_H1  * rho  * profile /(mass*PROTONMASS*bb))

    return tau

class MetalLines:
    """Generate metal line spectra from simulation snapshot"""
    def __init__(self,los_table, snapshot, cloudy_dir="/home/spb/codes/ArepoCoolingTables/tmp_spb/", nbins = 1024):
        #Get los table from group list
        f = h5py.File(snapshot)
        self.box = f["Header"].attrs["BoxSize"]
        self.hubble = f["Header"].attrs["HubbleParam"]
        self.atime = f["Header"].attrs["Time"]
        self.redshift = f["Header"].attrs["Redshift"]
        Omega0 = f["Header"].attrs["Omega0"]
        OmegaLambda = f["Header"].attrs["OmegaLambda"]
        self.Hz = 100.0*self.hubble * np.sqrt(Omega0/self.atime**3 + OmegaLambda)
        self.nbins = nbins
        self.xbins = np.arange(0,self.nbins)*self.box/self.nbins
        self.species = ['He', 'C', 'N', 'O', 'Ne', 'Mg', 'Si', 'Fe']
        self.NumLos = np.size(los_table.axis)
        #Line data
        self.lines = line_data.LineData()
        #generate metal and hydrogen spectral densities
        #Indexing is: rho_metals [ NSPECTRA, NBIN ]
        (self.rho_H, self.metals) = SPH_Interpolate_metals(f["PartType0"], los_table, nbins, self.box)
        #rescale H density
        self.rho_H = rescale_units_rho_H(self.rho_H, self.hubble, self.atime)
        #Rescale metals
        for (key, value) in self.metals:
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
        ion_density = np.array(species.rho)
        #Compute tau for this metal ion
        tau_metal=np.empty(np.shape(species.rho))
        for n in np.arange(0,self.NumLos):
            #For the density parameter use the hydrogen density at this pixel
            #For metallicity pass the metallicity of this species at this bin (rho_Z/ rho_H) and it will be converted to cloudy format
            for i in np.arange(0,self.nbins):
                ion_density[n,i] *= self.cloudy.ion(elem, ion, self.redshift, ion_density[n,i]/self.rho_H[n,i], self.rho_H[n,i])
            tau_metal[n] = compute_absorption(self.xbins, ion_density[n], species.vel[n], species.temp[n],line,self.Hz,self.hubble, self.box, self.atime,mass)
        return tau_metal


