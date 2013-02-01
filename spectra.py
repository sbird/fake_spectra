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

def SPH_Interpolate(data, los_table, nbins, box):
    """Interpolate particles to lines of sight, calculating density, temperature and velocity
    of various species along the line of sight.

    This is a wrapper which calls the C function.
    Arguments:
    	data - HDF5 dataset from snapshot. Use f["PartType0"]
	    los_table - table of los positions. should have member arrays x, y, z and axis.
	    nbins - number of bins in each spectrum
	    box - box size
        hub - hubble constant (eg 0.71)
        atime - 1/(1+z)

    Returns:
        (HI, metals)
        where HI is a Species class containing the neutral hydrogen
        and metals is a list of NMETALS Species classes, each containing a metal species
        If the GFM_Metals array is not found, will just return HI
    """
    pos = np.array(data["Coordinates"],dtype=np.float32)
    vel = np.array(data["Velocities"],dtype=np.float32)
    mass = np.array(data["Masses"],dtype=np.float32)
    u = np.array(data["InternalEnergy"],dtype=np.float32)
    nh0 = np.array(data["NeutralHydrogenAbundance"],dtype=np.float32)
    ne = np.array(data["ElectronAbundance"],dtype=np.float32)
    hh = np.array(hsml.get_smooth_length(data),dtype=np.float32)
    xx=np.array(los_table.xx, dtype=np.float32)
    yy=np.array(los_table.yy, dtype=np.float32)
    zz=np.array(los_table.zz, dtype=np.float32)
    axis=np.array(los_table.axis, dtype=np.int32)
    try:
        metal_in = np.array(data["GFM_Metals"],dtype=np.float32)[:,2:]
        #Deal with floating point roundoff - metal_in will sometimes be negative
        metal_in[np.where(np.abs(metal_in) < 1e-10)] = 0
        (rho_HI, vel_HI, temp_HI, rho_metal, vel_metal, temp_metal) =  _SPH_Interpolate(nbins, box, pos, vel, mass, u, nh0, ne, metal_in, hh, axis, xx, yy, zz)
        HI = Species(rho_HI, vel_HI, temp_HI)
        metals = [Species(rho_metal[:,:,0],vel_metal[:,:,0],temp_metal[:,:,0]),]
        for mm in np.arange(1,np.shape(metals)[1]):
            metals.append( Species(rho_metal[:,:,mm], vel_metal[:,:,mm], temp_metal[:,:,mm]))
        return (HI, metals)
    except IOError:
        metals = np.array()
        (rho_HI, vel_HI, temp_HI) =  _SPH_Interpolate(nbins, box, pos, vel, mass, u, nh0, ne, metals, hh, axis, xx, yy, zz)
        HI = Species(rho_HI, vel_HI, temp_HI)
        return HI


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
            Needed for the conversion between comoving kpc/h to physical m"""
        # Conversion factors from internal units
        rscale = (KPC*atime)/h100    # convert length to m
        vscale = np.sqrt(atime)        #convert velocity to kms^-1
        mscale = (1.0e10*SOLAR_MASS)/h100   # convert mass to kg
        escale = 1.0e6           # convert energy/unit mass to J kg^-1
        # convert (with mu) T to K
        tscale = ((GAMMA-1.0) * mass * PROTONMASS * escale ) / BOLTZMANN
        # Rescale density, vel and temp
        #unweight vel and temp by HI density
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
        self.species = ['C','N','O','Ne','Mg','Si','Fe']
        self.NumLos = np.size(los_table.axis)
        #Line data
        self.lines = line_data.LineData()
        #generate metal and hydrogen spectral densities
        #Indexing is: rho_metals [ NSPECTRA, NBIN ]
        (self.HI, self.metals) = SPH_Interpolate(f["PartType0"], los_table, nbins, self.box)
        masses = self.lines.masses.values()
        self.HI.rescale_units(self.hubble, self.atime, masses[0])
        for mm in np.arange(0, np.size(self.metals)):
            self.metals[mm].rescale_units(self.hubble, self.atime, masses[mm])
        #Compute tau for HI
        self.tau_HI=np.empty(np.shape(self.HI.rho))
        for n in np.arange(0,self.NumLos):
            self.tau_HI[n] = compute_absorption(self.xbins, self.HI.rho[n], self.HI.vel[n], self.HI.temp[n], self.lines.get_line('H',1),self.Hz,self.hubble, self.box, self.atime,self.lines.get_mass('H'))

        #Generate cloudy tables
        self.cloudy = convert_cloudy.CloudyTable(cloudy_dir)


    def get_lines(self, species, ion):
        """Get the optical depth for a particular species out of self.species:
           (C, N, O, Ne, Mg, Si, Fe)
           and some ion number
           NOTE: May wish to special-case SiIII at some point
        """
        spec_ind = self.species.index(species)
        metal_density = self.metals[spec_ind].rho
        #Use the total metallicity from summing metal species, not from the GFM_Metallicity
        #variable as the difference is small (~ 4%)
        tot_met = np.sum(metal_density,axis=0)
        ion = self.cloudy.interpolate(self.redshift, tot_met, metal_density, spec_ind, ion)
        #Generate density of this ion: cloudy densities are in log_10
        ion_density = 10**ion * metal_density
        #Compute tau for this metal ion
        tau_metal=np.empty(np.shape(self.metals[spec_ind].rho))
        for n in np.arange(0,self.NumLos):
            tau_metal[n] = compute_absorption(self.xbins, ion_density[n], self.metals[spec_ind].vel[n], self.metals[spec_ind].temp[n], self.lines.get_line(species,1),self.Hz,self.hubble, self.box, self.atime,self.lines.get_mass(species))
        return tau_metal


