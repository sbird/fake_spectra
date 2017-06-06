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
sum(GFM_Metals[:2]) +  GFM_Metallicity ~ 1

Also note that there is some instability at very low metallicities - the code will often return +-1e-20.
"""

from __future__ import print_function
import os
import os.path as path
import shutil
import numpy as np
import h5py

from . import abstractsnapshot as absn
from . import gas_properties
from . import line_data
from . import unitsystem
from . import voigtfit
from . import spec_utils
from . import fluxstatistics as fstat
from ._spectra_priv import _Particle_Interpolate, _near_lines

from .cloudy_tables import convert_cloudy
def _get_cloudy_table(red, cdir=None):
    """Get the cloudy table if we didn't already"""
    #Generate cloudy tables
    if cdir != None:
        return convert_cloudy.CloudyTable(red, cdir)
    else:
        return convert_cloudy.CloudyTable(red)

#python2 compat
try:
    xrange(1)
except NameError:
    xrange = range

class Spectra(object):
    """Class to interpolate particle densities along a line of sight and calculate their absorption
        Arguments:
            num - Snapshot number
            base - Name of base directory for snapshot
            cofm - table of los positions, as [n, 3] shape array.
            axis - axis along which to put the sightline
            res (optional) - Spectra pixel resolution in km/s
            snr - If nonzero, add noise for the requested signal to noise for the spectra, when loading from disc.
            spec_res - Resolution of the spectrograph. The spectra will be convolved with a Gaussian of this FWHM on loading from disc.
            load_halo - Should I load the halo catalogue?
            units - UnitSystem instance.
            sf_neutral - bug fix for certain Gadget versions. See gas_properties.py
    """
    def __init__(self,num, base,cofm, axis, res=1., cdir=None, savefile="spectra.hdf5", savedir=None, reload_file=False, snr = 0., spec_res = 8,load_halo=False, units=None, sf_neutral=True,quiet=False):
        #Present for compatibility. Functionality moved to HaloAssignedSpectra
        _= load_halo
        self.num = num
        self.base = base
        #Create the unit system
        if units != None:
            self.units = units
        else:
            self.units = unitsystem.UnitSystem()
        #Empty dictionary to add results to
        self.tau_obs = {}
        self.tau = {}
        self.sfr = {}
        self.vel_widths = {}
        self.absorber_width = {}
        self.colden = {}
        self.velocity = {}
        self.temp = {}
        #A cache of the indices of particles near sightlines.
        self.part_ind = {}
        #This variable should be set to true once the sightlines are fixed, and the cache can be used.
        self.cofm_final = False
        self.num_important = {}
        self.discarded=0
        self.npart=0
        #If greater than zero, will add noise to spectra when they are loaded.
        self.snr = snr
        self.spec_res = spec_res
        self.cdir = cdir
        #Minimum length of spectra within which to look at metal absorption (in km/s)
        self.minwidth = 500.
        try:
            self.snapshot_set = absn.AbstractSnapshotFactory(num, base)
        except IOError:
            pass
        if savedir is None:
            #Use snapdir if exists, otherwise use SPEC_
            savedir = path.join(base,"snapdir_"+str(num).rjust(3,'0'))
            #Make sure savedir exists.
            if not path.exists(savedir):
                savedir = path.join(base,"SPECTRA_"+str(num).rjust(3,'0'))
        self.savefile = path.join(savedir,savefile)
        #Snapshot data
        if reload_file:
            if not quiet:
                print("Reloading from snapshot (will save to: ",self.savefile," )")
            #Make sure the obvious syntax for a single sightline works
            if np.shape(cofm) == (3,):
                cofm = np.array([cofm,])
            self.cofm = cofm.astype(np.float64)
            if np.shape(axis) == ():
                axis = np.array([axis])
            self.axis = axis.astype(np.int32)
            if cofm is None or axis is None:
                raise RuntimeError("None was passed for cofm or axis. If you are trying to load from a savefile, use reload_file=False.")
            try:
                self.npart=self.snapshot_set.get_npart()
                #If we got here without a snapshot_set, we really have an IOError
            except AttributeError:
                raise IOError("Unable to load snapshot ",num, base)
            self.box = self.snapshot_set.get_header_attr("BoxSize")
            self.atime = self.snapshot_set.get_header_attr("Time")
            self.red = 1/self.atime - 1.
            self.hubble = self.snapshot_set.get_header_attr("HubbleParam")
            self.OmegaM = self.snapshot_set.get_header_attr("Omega0")
            self.OmegaLambda = self.snapshot_set.get_header_attr("OmegaLambda")
            #Calculate omega_baryon (approximately only for HDF5)
            self.omegab = self.snapshot_set.get_omega_baryon()
            #Get the unit system.
            try:
                self.units = self.snapshot_set.get_units()
            except KeyError:
                if not quiet:
                    print('No units found. Using kpc/kms/10^10Msun by default')
        else:
            self.load_savefile(self.savefile)
        # Conversion factors from internal units
        self.rscale = (self.units.UnitLength_in_cm*self.atime)/self.hubble
        # Convert length to physical cm
        # Calculate the length scales to be used in the box: Hz in km/s/Mpc
        Hz = 100.0*self.hubble * np.sqrt(self.OmegaM/self.atime**3 + self.OmegaLambda)
        #Convert comoving internal units to physical km/s.
        #Numerical constant is 1 Mpc in cm.
        self.velfac = self.rscale * Hz / 3.085678e24
        self.vmax = self.box * self.velfac # box size (physical kms^-1)
        self.NumLos = np.size(self.axis)
        try:
            # velocity bin size (kms^-1)
            self.dvbin = self.vmax / (1.*self.nbins)
        except AttributeError:
            #This will occur if we are not loading from a savefile
            self.dvbin = res # velocity bin size (kms^-1)
            #Number of bins to achieve the required resolution
            self.nbins = int(self.vmax / self.dvbin)
        #Species we can use: Z is total metallicity
        self.species = ['H', 'He', 'C', 'N', 'O', 'Ne', 'Mg', 'Si', 'Fe', 'Z']
        #Solar abundances from Asplund 2009 / Grevasse 2010 (which is used in Cloudy 13, Hazy Table 7.4).
        self.solar = {"H":1, "He":0.0851, "C":2.69e-4,"N":6.76e-5,"O":4.9e-4,"Ne":8.51e-5,"Mg":3.98e-5,"Si":3.24e-5,"Fe":3.16e-5}
        # Total solar metallicity is from Asplund 2009 0909.0948
        # Note the solar metallicity is the mass fraction of metals
        # divided by the mass fraction of hydrogen
        self.solarz = 0.0134/0.7381
        #Line data
        self.lines = line_data.LineData()
        #Load the class for computing gas properties such as temperature from the raw simulation.
        try:
            self.gasprop=gas_properties.GasProperties(redshift = self.red, absnap=self.snapshot_set, hubble=self.hubble, units=self.units, sf_neutral=sf_neutral)
        except AttributeError:
            #Occurs if we didn't load a snapshot
            pass
        if not quiet:
            print(self.NumLos, " sightlines. resolution: ", self.dvbin, " z=", self.red)

    def save_file(self):
        """
        Saves spectra to a file, because they are slow to generate.
        File is by default to be $snap_dir/snapdir_$snapnum/spectra.hdf5.
        """
        #We should make sure we have loaded all lazy-loaded things first.
        self._load_all_multihash(self.tau_obs, "tau_obs")
        self._load_all_multihash(self.tau, "tau")
        self._load_all_multihash(self.colden, "colden")
        try:
            self._load_all_multihash(self.velocity, "velocity")
        except IOError:
            pass
        #Make sure the directory exists
        if not path.exists(path.dirname(self.savefile)):
            os.mkdir(path.dirname(self.savefile))
        #Make a backup.
        if path.exists(self.savefile):
            shutil.move(self.savefile,self.savefile+".backup")
        try:
            f=h5py.File(self.savefile,'w')
        except IOError:
            raise IOError("Could not open ",self.savefile," for writing")
        self._save_file(f)

    def _save_file(self, f):
        """Saves to an open file handle, so it can be called by child classes which may want to save extra data."""
        grp = f.create_group("Header")
        grp.attrs["redshift"]=self.red
        grp.attrs["nbins"]=self.nbins
        grp.attrs["hubble"]=self.hubble
        grp.attrs["box"]=self.box
        grp.attrs["omegam"]=self.OmegaM
        grp.attrs["omegab"]=self.omegab
        grp.attrs["omegal"]=self.OmegaLambda
        grp.attrs["discarded"]=self.discarded
        grp.attrs["npart"]=self.npart
        grp = f.create_group("spectra")
        grp["cofm"]=self.cofm
        grp["axis"]=self.axis
        #Observer tau is the strongest unsaturated line
        grp_grid = f.create_group("tau_obs")
        self._save_multihash(self.tau_obs, grp_grid)
        #Optical depth in specific lines
        grp_grid = f.create_group("tau")
        self._save_multihash(self.tau, grp_grid)
        #Column density
        grp_grid = f.create_group("colden")
        self._save_multihash(self.colden, grp_grid)
        #Velocity
        grp_grid = f.create_group("velocity")
        self._save_multihash(self.velocity, grp_grid)
        #Temperature
        grp_grid = f.create_group("temperature")
        self._save_multihash(self.temp, grp_grid)
        #Number of particles important for each spectrum
        grp_grid = f.create_group("num_important")
        self._save_multihash(self.num_important, grp_grid)
        f.close()

    def _load_all_multihash(self,array, array_name):
        """Do all allowed lazy-loading for an array.
        """
        for key in list(array.keys()):
            self._really_load_array(key, array, array_name)

    def _save_multihash(self,save_array, grp):
        """Save an array using a tuple key, like save_array[(elem, ion, line)]
        to a hierarchy of hdf groups below grp"""
        for (key, value) in save_array.items():
            #Create directory hierarchy recursively
            gg = grp
            for ii in range(np.size(key)-1):
                try:
                    gg = gg[str(key[ii])]
                except KeyError:
                    gg.create_group(str(key[ii]))
                    gg = gg[str(key[ii])]
            #Delete old dataset if present
            try:
                del gg[str(key[-1])]
            except KeyError:
                pass
            #Save the dataset
            gg.create_dataset(str(key[-1]),data=value)

    def _really_load_array(self, key, array, array_name):
        """Replace a lazy-loaded array with the real one from disc"""
        #First check it was not already loaded
        if np.size(array[key]) > 1:
            return
        #If not, load it.
        f=h5py.File(self.savefile,'r')
        if np.size(key) == 2:
            array[key] = np.array(f[array_name][str(key[0])][str(key[1])])
        elif np.size(key) == 3:
            array[key] = np.array(f[array_name][str(key[0])][str(key[1])][str(key[2])])
        else:
            raise ValueError("Not supported")
        f.close()

    def add_noise(self, snr, tau, seed):
        """Compute a Gaussian noise vector from the flux variance and the SNR, as computed from optical depth"""
        flux = np.exp(-tau)
        if np.size(np.shape(flux)) == 1:
            lines = 1
        else:
            lines = np.shape(flux)[0]
        #This is to get around the type rules.
        if lines == 1:
            #This ensures that we always get the same noise for the same spectrum
            np.random.seed(seed)
            flux += np.random.normal(0, 1./snr, self.nbins)
        else:
            for ii in xrange(lines):
                np.random.seed(ii)
                flux[ii]+=np.random.normal(0,1./snr, self.nbins)
        #Make sure we don't have negative flux
        ind = np.where(flux > 0)
        tau[ind] = -np.log(flux[ind])
#         assert np.all(np.logical_not(np.isnan(tau)))
        return tau

    def load_savefile(self,savefile=None):
        """Load data from a file"""
        #Name of savefile
        try:
            f=h5py.File(savefile,'r')
        except IOError:
            raise IOError("Could not read saved data from: "+savefile+". If the file does not exist, try using reload_file=True")
        grid_file=f["Header"]
        self.red=grid_file.attrs["redshift"]
        self.atime = 1./(1+self.red)
        self.OmegaM=grid_file.attrs["omegam"]
        self.nbins=grid_file.attrs["nbins"]
        self.omegab=grid_file.attrs["omegab"]
        self.OmegaLambda=grid_file.attrs["omegal"]
        self.hubble=grid_file.attrs["hubble"]
        self.npart = np.array(grid_file.attrs["npart"])
        self.box=grid_file.attrs["box"]
        self.discarded=grid_file.attrs["discarded"]
        grp = f["colden"]
        for elem in grp.keys():
            for ion in grp[elem].keys():
                self.colden[(elem, int(ion))] = np.array([0])
        grp = f["tau_obs"]
        for elem in grp.keys():
            for ion in grp[elem].keys():
                self.tau_obs[(elem, int(ion))] = np.array([0])
        grp = f["tau"]
        for elem in grp.keys():
            for ion in grp[elem].keys():
                for line in grp[elem][ion].keys():
                    self.tau[(elem, int(ion),int(float(line)))] = np.array([0])
        try:
            grp = f["velocity"]
            for elem in grp.keys():
                for ion in grp[elem].keys():
                    self.velocity[(elem, int(ion))] = np.array([0])
        except KeyError:
            pass
        try:
            grp = f["temperature"]
            for elem in grp.keys():
                for ion in grp[elem].keys():
                    self.temp[(elem, int(ion))] = np.array([0])
        except KeyError:
            pass
        grp = f["num_important"]
        for elem in grp.keys():
            for ion in grp[elem].keys():
                self.num_important[(elem, int(ion))] = np.array(grp[elem][ion])
        grp=f["spectra"]
        self.cofm = np.array(grp["cofm"])
        self.axis = np.array(grp["axis"])
        f.close()

    def _interpolate_single_file(self,nsegment, elem, ion, ll, get_tau):
        """Read arrays and perform interpolation for a single file"""
        (pos, vel, elem_den, temp, hh, amumass) = self._read_particle_data(nsegment, elem, ion,get_tau)
        if amumass is False:
            return np.zeros([np.shape(self.cofm)[0],self.nbins],dtype=np.float32)
        if get_tau:
            #Allow us to compute absorption profiles assuming all
            #of one element is in the absorbing state
            if ion == -1:
                #Slight hack: try every ion in turn until we find the one with this line
                for ii in range(8):
                    try:
                        line = self.lines[(elem, ii)][ll]
                    except KeyError:
                        continue
                    break
            else:
                line = self.lines[(elem,ion)][ll]
        else:
            #If get_tau is false, we don't use the line data
            #(only the density), so we can set this as we like.
            #Setting something makes writing the interfaces easier,
            #because C doesn't have default arguments.
            line = self.lines[("H",1)][1215]
        return self._do_interpolation_work(pos, vel, elem_den, temp, hh, amumass, line, get_tau)

    def _read_particle_data(self,fn, elem, ion, get_tau):
        """Read the particle data for a single interpolation"""
        pos = self.snapshot_set.get_data(0,"Position",segment = fn).astype(np.float32)
        hh = self.snapshot_set.get_smooth_length(0,segment=fn).astype(np.float32)

        #Find particles we care about
        if self.cofm_final:
            try:
                ind = self.part_ind[fn]
            except KeyError:
                ind = self.particles_near_lines(pos, hh,self.axis,self.cofm)
                self.part_ind[fn] = ind
        else:
            ind = self.particles_near_lines(pos, hh,self.axis,self.cofm)
        #Do nothing if there aren't any, and return a suitably shaped zero array
        if np.size(ind) == 0:
            return (False, False, False, False,False,False)
        pos = pos[ind,:]
        hh = hh[ind]
        #Get the rest of the arrays: reducing them each time to have a smaller memory footprint
        vel = np.zeros(1,dtype=np.float32)
        temp = np.zeros(1,dtype=np.float32)
        if get_tau:
            vel = self.snapshot_set.get_data(0,"Velocity",segment = fn).astype(np.float32)
            vel = vel[ind,:]
        #gas density amu / cm^3
        den=self.gasprop.get_code_rhoH(0, segment=fn).astype(np.float32)
        # Get mass of atomic species
        if elem != "Z":
            amumass = self.lines.get_mass(elem)
        else:
            amumass = 1
        den = den[ind]
        #Only need temp for ionic density, and tau later
        if get_tau or (ion != -1 and elem != 'H'):
            temp = self.snapshot_set.get_temp(0,segment=fn, units=self.units).astype(np.float32)
            temp = temp[ind]
        #Find the mass fraction in this ion
        #Get the mass fraction in this species: elem_den is now density in ionic species in amu/cm^2 kpc/h
        #(these weird units are chosen to be correct when multiplied by the smoothing length)
        elem_den = (den*self.rscale)*self.get_mass_frac(elem,fn,ind)
        #Special case H1:
        if elem == 'H' and ion == 1:
            # Neutral hydrogen mass frac
            elem_den *= (self.gasprop.get_reproc_HI(0, segment=fn)[ind]).astype(np.float32)
        elif ion != -1:
            #Cloudy density in physical H atoms / cm^3
            ind2 = self._filter_particles(elem_den, pos, vel, den)
            if np.size(ind2) == 0:
                return (False, False, False, False,False,False)
            #Shrink arrays: we don't want to interpolate particles
            #with no mass in them
            temp = temp[ind2]
            pos = pos[ind2]
            hh = hh[ind2]
            if get_tau:
                vel = vel[ind2]
            elem_den = elem_den[ind2] * self._get_elem_den(elem, ion, den[ind2], temp, ind, ind2)
            del ind2
        #Get rid of ind so we have some memory for the interpolator
        del den
        #Put density into number density of particles, from amu
        elem_den/=amumass
        #Do interpolation.
        return (pos, vel, elem_den, temp, hh, amumass)

    def _filter_particles(self, elem_den, pos, velocity, den):
        """Get a filtered list of particles to add to the sightlines"""
        _ = (pos,velocity, den)
        ind2 = np.where(elem_den > 0)
        return ind2

    def _get_elem_den(self, elem, ion, den, temp, ind, ind2):
        """Get the density in an elemental species. Broken out so it can be over-ridden by child classes."""
        #Shut up a pylint warning
        _ = (ind, ind2)
        #Load a cloudy table if not done already
        try:
            self.cloudy_table
        except AttributeError:
            self.cloudy_table = _get_cloudy_table(self.red, self.cdir)
        #Make sure temperature doesn't overflow the cloudy table
        #High temperatures are unlikely to be in ionisation equilibrium anyway.
        #Low temperatures can be neglected because we don't follow cooling processes that far anyway.
        tmplimits = self.cloudy_table.get_temp_bounds()
        if np.max(temp) > tmplimits[1] or np.min(temp) < tmplimits[0]:
            temp2 = np.array(temp)
            temp2[np.where(temp2 > tmplimits[1])] = tmplimits[1]
            temp2[np.where(temp2 < tmplimits[0])] = tmplimits[0]
        else:
            temp2 = temp
        #Ditto density: high densities are not followed correctly and low densities contain no stuff.
        denslimits = self.cloudy_table.get_dens_bounds()
        if np.max(den) > denslimits[1] or np.min(den) < denslimits[0]:
            den2 = np.array(den)
            den2[np.where(den2 > denslimits[1])] = denslimits[1]
            den2[np.where(den2 < denslimits[0])] = denslimits[0]
        else:
            den2 = den
        return np.float32(self.cloudy_table.ion(elem, ion, den2, temp2))

    def _do_interpolation_work(self,pos, vel, elem_den, temp, hh, amumass, line, get_tau):
        """Run the interpolation on some pre-determined arrays, spat out by _read_particle_data"""
        #Factor of 10^-8 converts line width (lambda_X) from Angstrom to cm
        return _Particle_Interpolate(get_tau*1, self.nbins, self.snapshot_set.get_kernel(), self.box, self.velfac, self.atime, line.lambda_X*1e-8, line.gamma_X, line.fosc_X, amumass, pos, vel, elem_den, temp, hh, self.axis, self.cofm)

    def particles_near_lines(self, pos, hh,axis=None, cofm=None):
        """Filter a particle list, returning an index list of those near sightlines"""
        if axis is None:
            axis = self.axis
        if cofm is None:
            cofm = self.cofm
        #Axis is 1-indexed between 1 and 3. 1 is x axis.
        assert np.min(axis) > 0
        assert np.max(axis) <4
        ind = _near_lines(self.box, pos, hh, axis, cofm)
        return ind

    def get_mass_frac(self,elem, fn, ind):
        """Get the mass fraction of a given species from a snapshot.
        Arguments:
            elem = name of element
            data = pointer to hdf5 array containing baryons
            ind = index of particles we care about
        Returns mass_frac - mass fraction of this ion
        """
        if elem == "Z":
            mass_frac = self.snapshot_set.get_data(0,"Metallicity",segment = fn).astype(np.float32)
        else:
            nelem = self.species.index(elem)
            #Get metallicity of this metal species
            try:
                mass_frac = (self.snapshot_set.get_data(0,"GFM_Metals",segment = fn).astype(np.float32))[:,nelem]
            except KeyError:
                #If GFM_Metals is not defined, fall back to primordial abundances
                metal_abund = np.array([0.76, 0.24],dtype=np.float32)
                nvalues = self.snapshot_set.get_blocklen(0,"Density",segment = fn)
                mass_frac = metal_abund[nelem]*np.ones(nvalues, dtype=np.float32)
        mass_frac = mass_frac[ind]
        #Deal with floating point roundoff - mass_frac will sometimes be negative
        mass_frac[np.where(mass_frac <= 0)] = 0
        assert mass_frac.dtype == np.float32
        return mass_frac

    def replace_not_DLA(self, ndla, thresh=10**20.3, elem="H", ion=1):
        """
        Replace those sightlines which do not contain sightlines above a given column density with new sightlines, until all sightlines are above the column density.
        Keep track of the number discarded in self.discarded.
        Must implement get_cofm for this to work
        """
        #Declare variables
        found = 0
        wanted = ndla
        cofm_DLA = np.empty_like(self.cofm)[:ndla, :]
        #Filter
        #Note: line does nothing
        col_den = self.compute_spectra(elem,ion,1215,False)
        ind = self.filter_DLA(col_den, thresh)
        H1_DLA = np.empty_like(col_den)
        #Update saves
        top = np.min([wanted, found+np.size(ind)])
        cofm_DLA[found:top] = self.cofm[ind][:top,:]
        H1_DLA[found:top] = col_den[ind][:top,:]
        found += np.size(ind)
        self.discarded = self.NumLos-np.size(ind)
        print("Discarded: ",self.discarded)
        while found < wanted:
            #Get a bunch of new spectra
            self.cofm = self.get_cofm()
            col_den = self.compute_spectra(elem,ion,1215,False)
            ind = self.filter_DLA(col_den, thresh)
            #Update saves
            top = np.min([wanted, found+np.size(ind)])
            cofm_DLA[found:top] = self.cofm[ind][:top-found,:]
            H1_DLA[found:top] = col_den[ind][:top-found,:]
            found += np.size(ind)
            self.discarded += self.NumLos-np.size(ind)
            print("Discarded: ",self.discarded)
        #Correct proportions in case we find slightly more than we need
        self.discarded = int(self.discarded*1.*wanted/1./found)
        #Copy back
        self.cofm=cofm_DLA
        self.axis = self.axis[:ndla]
        self.colden[("H",1)]=H1_DLA
        #Finalise the cofm array
        self.cofm_final = True
        self.NumLos = ndla

    def get_cofm(self, num = None):
        """Find a bunch more sightlines: should be overridden by child classes"""
        raise NotImplementedError

    def filter_DLA(self, col_den, thresh=10**20.3):
        """Find sightlines with a DLA"""
        #DLAs are huge objects in redshift space (several 10s of A wide), so we want to
        #sum the column densities over the entire spectrum.
        cdsum = np.sum(col_den, axis=1)
        if np.size(thresh) > 1:
            ind = np.where(np.logical_and(cdsum > thresh[0], cdsum < thresh[1]))
        else:
            ind = np.where(cdsum > thresh)
        return ind

    def get_metallicity(self, width=0.):
        """Return the metallicity, as M/H.
        If width > 0, computes M/H +- width km/s from the maximum H peak."""
        MM = self.get_density("Z",-1)
        HH = self.get_density("H",-1)
        if width > 0:
            (roll, hhr) = spec_utils.get_rolled_spectra(HH)
            mmr = np.array([np.roll(MMr, rr) for (MMr,rr) in zip(MM,roll)])
            imax = int(np.shape(MM)[1]/2)
            mms = np.array([np.sum(mmrr[imax-width/self.dvbin:imax+width/self.dvbin]) for mmrr in mmr])
            hhs = np.array([np.sum(hhrr[imax-width/self.dvbin:imax+width/self.dvbin]) for hhrr in hhr])
        else:
            mms = np.sum(MM, axis=1)
            hhs = np.sum(HH, axis=1)
        return mms/hhs/self.solarz
        #Use only DLA regions: tricky to preserve shape
        #ma_HH = np.ma.masked_where(HH < thresh, MM/HH)
        #data = np.array([np.mean(ma_HH, axis=1)])
        #return data/solar

    def get_ion_metallicity(self, species,ion):
        """Get the metallicity derived from an ionic species"""
        MM = self.get_density(species,ion)
        HH = self.get_density("H",1)
        mms = np.sum(MM, axis=1)
        hhs = np.sum(HH, axis=1)
        return mms/hhs/self.solar[species]

    def compute_spectra(self,elem, ion, ll, get_tau):
        """
        Computes the absorption spectrum (tau (u) ) from a binned set of interpolated
        densities, velocities and temperatures.

        Optical depth is given by:
        tau (u) = sigma_X c / H(z) int_infty^infty n_x(x) V( u - x - v_pec, b(x) ) dx
        where V is the Voigt profile, b(x)^2 = 2k_B T /m_x c^2 is the velocity dispersion.
        and v_pec is the peculiar velocity.
        sigma_X is the cross-section for this transition.
        """
        #Get array sizes
        nsegments = self.snapshot_set.get_n_segments()
        result =  self._interpolate_single_file(0, elem, ion, ll, get_tau)
        #Do remaining files
        if nsegments == 1:
            return result
        for nn in xrange(1,nsegments):
            tresult =  self._interpolate_single_file(nn, elem, ion, ll, get_tau)
            #Add new file
            result += tresult
            del tresult
        return result

    def equivalent_width(self, elem, ion, line):
        """Calculate the equivalent width of a line in Angstroms"""
        tau = self.get_tau(elem, ion, line)
        #1 bin in wavelength: δλ =  λ . v / c
        #λ here is the rest wavelength of the line.
        #speed of light in km /s
        light = self.units.light / 1e5
        #lambda in Angstroms, dvbin in km/s,
        #so dl is in Angstrom
        dl = self.dvbin / light * line
        eq_width = np.trapz(-np.expm1(-tau),dx=dl, axis=1)
        #Don't need to divide by 1+z as lambda_X is already rest wavelength
        return eq_width

    def eq_width_hist(self, elem, ion, line, dv=0.05):
        """
        Compute a histogram of the equivalent width distribution of our spectra.

        Returns:
            (v, f_table) - v (binned in log) and corresponding f(N)
        """
        print("For ",line," Angstrom")
        eq_width = self.equivalent_width(elem, ion, line)
        v_table = np.arange(np.log10(np.min(eq_width)), np.log10(np.max(eq_width)), dv)
        vbin = np.array([(v_table[i]+v_table[i+1])/2. for i in range(0,np.size(v_table)-1)])
        eqws = np.histogram(np.log10(eq_width),v_table, density=True)[0]
        return (vbin, eqws)

    def get_col_density(self, elem, ion, force_recompute=False):
        """get the column density in each pixel for a given species.
        In units of [metal] ions cm^-2."""
        try:
            if force_recompute:
                raise KeyError
            self._really_load_array((elem, ion), self.colden, "colden")
            return self.colden[(elem, ion)]
        except KeyError:
            colden = self.compute_spectra(elem, ion, 0,False)
            self.colden[(elem, ion)] = colden
            return colden

    def get_density(self, elem, ion, force_recompute=False):
        """Get the density in each pixel for a given species, in units of [metal] ions cm^-3."""
        colden = self.get_col_density(elem, ion, force_recompute)
        phys = self.dvbin/self.velfac*self.rscale
        return colden/phys

    def get_tau(self, elem, ion,line, number = -1, force_recompute=False, noise=True):
        """Get the optical depth in each pixel along the sightline for a given line."""
        try:
            if force_recompute:
                raise KeyError
            self._really_load_array((elem, ion, line), self.tau, "tau")
            tau = self.tau[(elem, ion,line)]
        except KeyError:
            tau = self.compute_spectra(elem, ion, line,True)
            self.tau[(elem, ion,line)] = tau
        if number >= 0:
            tau = tau[number,:]
        tau = spec_utils.res_corr(tau, self.dvbin, self.spec_res)
        if noise and self.snr > 0:
            tau = self.add_noise(self.snr, tau, number)
        return tau

    def _vel_single_file(self,fn, elem, ion):
        """Get the column density weighted interpolated velocity field for a single file"""
        (pos, vel, elem_den, temp, hh, amumass) = self._read_particle_data(fn, elem, ion,True)
        if amumass is False:
            return np.zeros([np.shape(self.cofm)[0],self.nbins,3],dtype=np.float32)
        else:
            line = self.lines[("H",1)][1215]
            vv =  np.empty([np.shape(self.cofm)[0],self.nbins,3],dtype=np.float32)
            phys = self.dvbin/self.velfac*self.rscale
            for ax in (0,1,2):
                weight = vel[:,ax]*np.sqrt(self.atime)
                vv[:,:,ax] = self._do_interpolation_work(pos, vel, elem_den*weight/phys, temp, hh, amumass, line, False)
            return vv

    def _get_mass_weight_quantity(self, func, elem, ion):
        """
        Helper function to get a mass weighted quantity, which reduces code duplication.
        func should be something like _vel_single_file above (for velocity)
        and have the signature func(file, elem, ion)
        """
        nsegments = self.snapshot_set.get_n_segments()
        result =  func(0, elem, ion)
        if nsegments > 1:
            #Do remaining files
            for nn in xrange(1, nsegments):
                tresult =  func(nn, elem, ion)
                #Add new file
                result += tresult
                del tresult
        den = self.get_density(elem, ion)
        den[np.where(den == 0.)] = 1
        try:
            result /= den
        except ValueError:
            #Broadcasting can't handle velocity as it is 3d/1d
            for ax in range(3):
                result[:,:,ax] /= den
        return result

    def get_velocity(self, elem, ion):
        """Get the column density weighted velocity in each pixel for a given species.
        """
        try:
            self._really_load_array((elem, ion), self.velocity, "velocity")
            return self.velocity[(elem, ion)]
        except KeyError:
            velocity =  self._get_mass_weight_quantity(self._vel_single_file, elem, ion)
            self.velocity[(elem, ion)] = velocity
            return velocity

    def _temp_single_file(self,fn, elem, ion):
        """Get the density weighted interpolated temperature field for a single file"""
        (pos, vel, elem_den, temp, hh, amumass) = self._read_particle_data(fn, elem, ion,True)
        if amumass is False:
            return np.zeros([np.shape(self.cofm)[0],self.nbins],dtype=np.float32)
        else:
            line = self.lines[("H",1)][1215]
            phys = self.dvbin/self.velfac*self.rscale
            temp = self._do_interpolation_work(pos, vel, elem_den*temp/phys, temp, hh, amumass, line, False)
            return temp

    def get_temp(self, elem, ion):
        """Get the density weighted temperature in each pixel for a given species.
        """
        try:
            self._really_load_array((elem, ion), self.temp, "temperature")
            return self.temp[(elem, ion)]
        except KeyError:
            temp =  self._get_mass_weight_quantity(self._temp_single_file, elem, ion)
            self.temp[(elem, ion)] = temp
            return temp

    def get_b_param_dist(self, elem="H", ion=1, line=1215):
        """Get the power law betweeen the 'minimum' b parameter and column density."""
        tau = self.get_tau(elem, ion, line)
        (b_0, gamm1) = voigtfit.get_b_param_dist(tau, self.dvbin, elem=elem, ion=ion, line=line)
        return (b_0, gamm1)

    def column_density_from_voigt(self, elem="H",ion=1,line=1215,dlogN=0.2,minN=13,maxN=23, close=50.,dX=True, nspectra=-1):
        """This computes the column density function using column densities from a Voigt profile fit.
        Concatenate objects closer than close km/s."""
        NHI_table = 10**np.arange(minN, maxN, dlogN)
        center = np.array([(NHI_table[i]+NHI_table[i+1])/2. for i in range(0,np.size(NHI_table)-1)])
        width =  np.array([NHI_table[i+1]-NHI_table[i] for i in range(0,np.size(NHI_table)-1)])
        #Number of lines
        tot_lines = self.NumLos+self.discarded
        #Absorption distance for each line
        if dX:
            dist=self.units.absorption_distance(self.box, self.red)
        else:
            dist=self.units.redshift_distance(self.box, self.red, self.OmegaM)
        tau = self.get_tau(elem, ion, line)
        if nspectra > 0:
            tau = tau[:nspectra, :]
            tot_lines *= nspectra/self.NumLos
        #Combine all lines closer than the close value into one.
        n_vals = voigtfit.get_voigt_systems(tau, self.dvbin, elem=elem,ion=ion, line=line,verbose=False, close=close)
        (tot_f_N, NHI_table) = np.histogram(n_vals,NHI_table)
        #The normalisation should be the total sightline distance.
        tot_f_N=tot_f_N/(width*dist*tot_lines)
        return (center, tot_f_N)

    def column_density_function(self,elem = "H", ion = 1, dlogN=0.2, minN=13, maxN=23., line=True, close=50.,dX=True):
        """
        This computes the absorber column density distribution function, which is the number
        of absorbers per sight line with column densities in the interval
        [N, N+dN] at the absorption distance X.
        The path can be either a whole sightline across the box, if line=True, or, if line=False,
        it can be a fixed spacing in the hubble flow.

        So we have f(N) = d n_abs/ dN dX
        and n_DLA(N) = number of absorbers per sightline in this column density bin.
                     1 sightline is defined to be one grid cell.
                     So this is (cells in this bins) / (no. of cells)
        ie, f(N) = n_abs / ΔN / ΔX
        Note f(N) has dimensions of cm^2, because N has units of cm^-2 and X is dimensionless.
        Note column density is number of *ions* per cm^2, not amu per cm^2.

        Parameters:
            dlogN - bin spacing
            minN - minimum log N
            maxN - maximum log N
            line - cddf for whole line or for each cell.
            close - amalgamate absorbers closer than X km/s
            dX - If true return dn/dN dX, if false return dn/dN dz
        Returns:
            (NHI, f_N_table) - N_HI (binned in log) and corresponding f(N)
        """
        NHI_table = 10**np.arange(minN, maxN, dlogN)
        center = np.array([(NHI_table[i]+NHI_table[i+1])/2. for i in range(0,np.size(NHI_table)-1)])
        width =  np.array([NHI_table[i+1]-NHI_table[i] for i in range(0,np.size(NHI_table)-1)])
        #Number of lines
        tot_lines = self.NumLos+self.discarded
        #Absorption distance for each line
        if dX:
            dist=self.units.absorption_distance(self.box, self.red)
        else:
            dist=self.units.redshift_distance(self.box, self.red, self.OmegaM)
        if line:
            rho = np.sum(self.get_col_density(elem, ion), axis=1)
        else:
            rho = self.get_col_density(elem, ion)
            cbins = np.max((int(np.round((close/self.dvbin))),1))
            rhob = np.array([np.sum(rho[:,cbins*i:cbins*(i+1)],axis=1) for i in xrange(int(np.shape(rho)[1]/cbins))]).T
            #Check that fp roundoff is not too severe: this can sometimes trigger for silly reasons
            #assert np.abs((np.sum(rhob) / np.sum(rho))-1) < 5e-2
            rho = rhob
        (tot_f_N, NHI_table) = np.histogram(rho,NHI_table)
        #The normalisation should be the total sightline distance.
        tot_f_N=tot_f_N/(width*dist*tot_lines)
        return (center, tot_f_N)

    def _rho_abs(self, thresh=10**20.3, upthresh=None, elem = "H", ion = 1):
        """Compute rho_abs, the sum of the mass in an absorber,
           divided by the volume of the spectra in g/cm^3 (comoving).
            Omega_DLA = m_p * avg. column density / (1+z)^2 / length of column
        """
        #Column density of ion in atoms cm^-2 (physical)
        col_den = np.sum(self.get_col_density(elem, ion), axis=1)
        if thresh > 0 or upthresh != None:
            HIden = np.sum(col_den[np.where((col_den > thresh)*(col_den < upthresh))])/np.size(col_den)
        else:
            HIden = np.mean(col_den)
        HIden *= np.size(col_den)/(np.size(col_den)+1.*self.discarded)
        #Avg. Column density in g cm^-2 (comoving)
        HIden = self.lines.get_mass(elem) * self.units.protonmass * HIden/(1+self.red)**2
        #Length of column (each cell) in comoving cm
        length = (self.box*self.units.UnitLength_in_cm/self.hubble)
        #Avg density in g/cm^3 (comoving)
        return HIden/length

    def rho_DLA(self, thresh=10**20.3):
        """Compute rho_DLA, the average density in DLAs. This is almost the same as the average density in HI.
           Units are 10^8 M_sun / Mpc^3 (comoving), like 0811.2003
        """
        #Avg density in g/cm^3 (comoving)
        rho_DLA = self._rho_abs(thresh)
        # 1e8 M_sun/Mpc^3 in g/cm^3
        conv = 0.01 * self.units.UnitMass_in_g / self.units.UnitLength_in_cm**3
        return rho_DLA / conv

    def omega_abs(self, thresh=10**20.3, upthresh=1e40, elem = "H", ion = 1):
        """Compute Omega_abs, the sum of the mass in a given absorber,
            divided by the volume of the spectra, divided by the critical density.
            Ω_abs = m_p * avg. column density / (1+z)^2 / length of column / rho_c
            Note: If we want the neutral gas density rather than the neutral hydrogen density, divide by 0.76,
            the hydrogen mass fraction.
        """
        #Avg density in g/cm^3 (comoving) divided by critical density in g/cm^3
        omega_DLA=self._rho_abs(thresh, upthresh, elem=elem, ion=ion)/self.units.rho_crit(self.hubble)
        return omega_DLA

    def omega_abs_cddf(self, thresh=10**20.3, upthresh=1e40, elem = "H", ion = 1):
        """Compute Omega_abs, the sum of the mass in a given absorber,
            divided by the volume of the spectra, divided by the critical density.
            Omega_abs = m_p * avg. column density / (1+z)^2 / length of column / rho_c
            This is computed by summing the column density function, rather than directly by summing
            columns. Should be the same as omega_abs.
        """
        dlogN = 0.2
        (bins, cddf) = self.column_density_function(elem, ion, dlogN, minN=np.log10(thresh), maxN=np.log10(upthresh))
        #Integrate cddf * N
        moment = cddf*bins
        #H0 in 1/s units
        h100=self.units.h100*self.hubble
        #The 1+z factor converts lightspeed to comoving
        omega_abs = self.lines.get_mass(elem)*self.units.protonmass/self.units.light*h100/self.units.rho_crit(self.hubble)*np.trapz(moment, bins)
        return omega_abs

    def line_density(self, thresh=10**20.3, upthresh=10**40, elem = "H", ion = 1):
        """Compute the line density, the total no. of DLA sightlines divided by the total number of sightlines, multiplied by d L / dX. This is dN/dX = l_DLA(z)
        """
        col_den = np.sum(self.get_col_density(elem, ion), axis=1)
        #Average fraction of pixels containing a DLA
        frac = 1.*np.size(col_den[np.where((col_den > thresh)*(col_den < upthresh))])/np.size(col_den)
        #Divide by abs. distance per sightline
        frac *= np.size(col_den)/(np.size(col_den)+1.*self.discarded)
        return frac/(self.units.absorption_distance(self.box, self.red))

    def line_density_eq_w(self, thresh=0.4, elem = "H", ion = 1, line=1216):
        """Compute the line density with an equivalent width threshold.
        This is the total no. of sightlines above the threshold divided by the total number of sightlines,
        multiplied by d L / dX. This is dN/dX = l(z)
        """
        eq_width = self.equivalent_width(elem, ion, line)
        #Average fraction of pixels containing a DLA
        frac = 1.*np.size(eq_width[np.where(eq_width > thresh)])
        frac /= (np.size(eq_width)+1.*self.discarded)
        #Divide by abs. distance per sightline
        return frac/self.units.absorption_distance(self.box, self.red)

    def get_spectra_proj_pos(self, cofm=None):
        """Get the position of the spectra projected to their origin"""
        if np.mean(self.axis) != self.axis[0] or  self.axis[0] != self.axis[-1]:
            raise ValueError("Not all spectra are along the same axis")
        if cofm is None:
            cofm = self.cofm
        axis = self.axis[0]
        if axis == 1:
            spos = cofm[:,1:]
        if axis == 2:
            spos = np.vstack([cofm[:,0],cofm[:,2]]).T
        if axis == 3:
            spos = cofm[:,:2]
        return spos

    def get_mean_flux(self, elem="H", ion=1, line=1215):
        """Get the mean flux along a set of sightlines"""
        tau = self.get_tau(elem, ion, line)
        return np.mean(np.exp(-tau))

    def get_flux_pdf(self, elem="H", ion=1, line=1215, nbins=20, mean_flux_desired=None):
        """Get the flux PDF, a histogram of the flux values."""
        tau = self.get_tau(elem, ion, line)
        return fstat.flux_pdf(tau, nbins=nbins, mean_flux_desired = mean_flux_desired)

    def get_flux_power_1D(self, elem="H",ion=1, line=1215, mean_flux_desired = None):
        """Get the power spectrum of (variations in) the flux along the line of sight.
        This is: P_F(k_F) = <d_F d_F>
                 d_F = e^-tau / mean(e^-tau) - 1
        If mean_flux_desired is set, the spectral optical depths will be rescaled
        to match the desired mean flux."""
        tau = self.get_tau(elem, ion, line)
        (kf, avg_flux_power) = fstat.flux_power(tau, self.vmax, mean_flux_desired=mean_flux_desired)
        return kf[1:],avg_flux_power[1:]
