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
    if cdir is None:
        return convert_cloudy.CloudyTable(red)
    return convert_cloudy.CloudyTable(red, cdir)

#python2 compat
try:
    xrange(1)
except NameError:
    xrange = range

class Spectra:
    """Class to interpolate particle densities along a line of sight and calculate their absorption
        Arguments:
            num - Snapshot number
            base - Name of base directory for snapshot
            cofm - table of los positions, as [n, 3] shape array. Can be None if reloading from a savefile.
            axis - axis along which to put the sightline. Can be None if reloading from a savefile.
        Optional arguments:
            res - Pixel width of the spectrum in km/s
            spec_res - Resolution of the simulated spectrograph in km/s. Note this is not the pixel width.
                       Spectra will be convolved with a Gaussian of this rms on loading from disc.
            cdir - Directory containing cloudy tables.
            savefile - name of file to save spectra to.
            savedir - Directory of file to save spectra to.
            reload_file - if true, ignore save spectra and regenerate from the snapshot.
            load_halo - Whether to load the halo catalogue
            units - If not None, UnitSystem instance which overrides the default units read from the simulation.
            sf_neutral - bug fix for certain Gadget versions. See gas_properties.py
            quiet - Whether to output debug messages
            load_snapshot - Whether to load the snapshot
            gasprop - Name (not instance!) of class to compute neutral fractions and temperatures.
                      It should inherit from gas_properties.GasProperties and provide get_reproc_HI
                      for neutral fractions and get_temp for temperatures.
                      Default is gas_properties.GasProperties which reads both of these from the particle output.
            gasprop_args - Dictionary of extra arguments to be fed to gasprop, if gasprop is not the default.
            kernel - Interpolation method to use. The default is to select the method based on the type of simulation:
                     this will use a Voronoi mesh build if we are Arepo and an SPH kernel for Gadget.
                     Other options are "voronoi" which forces the Voronoi kernel, "tophat" which forces a flat tophat
                     kernel (a good back up for large Arepo simulations) "quintic" for a quintic SPh kernel as used in modern SPH
                     and "cubic" or "sph" for an old-school cubic SPH kernel.
    """
    def __init__(self, num, base, cofm, axis, MPI=None, res=1., cdir=None, savefile="spectra.hdf5", savedir=None, reload_file=False, spec_res = 0,load_halo=False, units=None, sf_neutral=True,quiet=False, load_snapshot=True, gasprop=None, gasprop_args=None, kernel=None):

        #Present for compatibility. Functionality moved to HaloAssignedSpectra
        _= load_halo
        self.num = num
        self.base = base
        # Adpot the MPI communiactor if desired
        self.MPI = MPI
        if MPI is not None:
            self.comm = MPI.COMM_WORLD
            self.rank = self.comm.Get_rank()
            self.size = self.comm.Get_size()
        else :
            self.comm = None
            self.rank = 0
            self.size = 1
        #Create the unit system
        if units is not None:
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
        self.dens_weight_dens = {}
        #A cache of the indices of particles near sightlines.
        self.part_ind = {}
        #This variable should be set to true once the sightlines are fixed, and the cache can be used.
        self.cofm_final = False
        self.num_important = {}

        self.discarded=0
        self.npart=0

        self.spec_res = spec_res
        self.cdir = cdir
        #Minimum length of spectra within which to look at metal absorption (in km/s)
        self.minwidth = 500.
        #Stopping criterion for optical depth: if a particle is adding less
        #than this to a pixel, stop the integration.
        self.tautail = 1e-7
        try:
            if load_snapshot:
                self.snapshot_set = absn.AbstractSnapshotFactory(num, base, comm=self.comm)
                #Set up the kernel
                if kernel is None:
                    self.kernel_int = self.snapshot_set.get_kernel()
                elif kernel == "voronoi":
                    self.kernel_int = 2
                elif kernel == "tophat":
                    self.kernel_int = 0
                elif kernel == "quintic":
                    self.kernel_int = 3
                elif kernel == "cubic":
                    #Cubic SPH kernel
                    self.kernel_int = 1
                elif kernel == "sph":
                    self.kernel_int = 1
                else:
                    raise ValueError("Unrecognised kernel %d" % kernel)
        except IOError:
            pass
        if savedir is None:
            #Use snapdir if exists, otherwise use SPEC_
            savedir = path.join(base, "snapdir_"+str(num).rjust(3, '0'))
            #Make sure savedir exists.
            if not path.exists(savedir):
                savedir = path.join(base, "SPECTRA_"+str(num).rjust(3, '0'))
        self.savefile = path.join(savedir, savefile)
        #Snapshot data
        if reload_file:
            if not quiet:
                print("Reloading from snapshot (will save to: ", self.savefile, " )", flush=True)
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
                self.npart = self.snapshot_set.get_npart()
                #If we got here without a snapshot_set, we really have an IOError
            except AttributeError:
                raise IOError("Unable to load snapshot ", num, base)
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
        # specify pixel width and number of pixels in spectra
        if reload_file:
            # if reloading from snapshot, pixel width must have been defined
            if res is None:
                raise ValueError('pixel width (res) not provided')
            # nbins must be an integer
            self.nbins=int(self.vmax/res)
            # velocity bin size (kms^-1)
            self.dvbin = self.vmax / (1.*self.nbins)
        else:
            # self.nbins already set from file
            self.dvbin = self.vmax / (1.*self.nbins)
            # check if you asked for different pixel width
            if res is not None:
                assert np.isclose(self.dvbin,res,rtol=1e-2),'pixel width error'
        #Species we can use: Z is total metallicity
        self.species = ['H', 'He', 'C', 'N', 'O', 'Ne', 'Mg', 'Si', 'Fe', 'Z']
        #Solar abundances from Asplund 2009 / Grevasse 2010 (which is used in Cloudy 13, Hazy Table 7.4).
        self.solar = {"H":1, "He":0.0851, "C":2.69e-4, "N":6.76e-5, "O":4.9e-4, "Ne":8.51e-5, "Mg":3.98e-5, "Si":3.24e-5, "Fe":3.16e-5}
        # Total solar metallicity is from Asplund 2009 0909.0948
        # Note the solar metallicity is the mass fraction of metals
        # divided by the mass fraction of hydrogen
        self.solarz = 0.0134/0.7381
        #Line data
        self.lines = line_data.LineData()
        #Load the class for computing gas properties such as temperature from the raw simulation.
        if gasprop is None:
            gasprop = gas_properties.GasProperties
        try:
            gprop_args = {"redshift" : self.red, "absnap" : self.snapshot_set, "hubble": self.hubble, "units": self.units, "sf_neutral": sf_neutral}
            if gasprop_args is not None:
                gprop_args.update(gasprop_args)
            self.gasprop = gasprop(**gprop_args)
        except AttributeError:
            #Occurs if we didn't load a snapshot
            pass
        if not quiet:
            print(self.NumLos, " sightlines. resolution: ", self.dvbin, " z=", self.red)

    def save_file(self):
        """
        Saves spectra to a file, because they are slow to generate.
        File is by default to be $snap_dir/snapdir_$snapnum/spectra.hdf5.
        Rank = 0 saves the full spectra.
        """
        if self.rank == 0:
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
                shutil.move(self.savefile, self.savefile+".backup")
            try:
                f = h5py.File(self.savefile, 'w')
            except IOError:
                raise IOError("Could not open ", self.savefile, " for writing")
            self._save_file(f)

    def _save_file(self, f):
        """Saves to an open file handle, so it can be called by child classes which may want to save extra data."""
        grp = f.create_group("Header")
        grp.attrs["redshift"] = self.red
        grp.attrs["nbins"] = self.nbins
        grp.attrs["hubble"] = self.hubble
        grp.attrs["box"] = self.box
        grp.attrs["omegam"] = self.OmegaM
        grp.attrs["omegab"] = self.omegab
        grp.attrs["omegal"] = self.OmegaLambda
        grp.attrs["discarded"] = self.discarded
        grp.attrs["npart"] = self.npart
        grp = f.create_group("spectra")
        grp["cofm"] = self.cofm
        grp["axis"] = self.axis
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
        #Density weighted density for each spectrum:
        #see get_dens_weighted_density for a description of what this is.
        grp_grid = f.create_group("density_weight_density")
        self._save_multihash(self.dens_weight_dens, grp_grid)
        f.close()

    def _load_all_multihash(self, array, array_name):
        """Do all allowed lazy-loading for an array.
        """
        for key in list(array.keys()):
            self._really_load_array(key, array, array_name)

    def _save_multihash(self, save_array, grp):
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
            gg.create_dataset(str(key[-1]), data=value)

    def _really_load_array(self, key, array, array_name):
        """Replace a lazy-loaded array with the real one from disc"""
        #First check it was not already loaded
        if np.size(array[key]) > 1:
            return
        #If not, load it.
        f = h5py.File(self.savefile, 'r')
        if np.size(key) == 2:
            array[key] = np.array(f[array_name][str(key[0])][str(key[1])])
        elif np.size(key) == 3:
            array[key] = np.array(f[array_name][str(key[0])][str(key[1])][str(key[2])])
        else:
            raise ValueError("Not supported")
        f.close()

    def add_noise(self, snr, flux, spec_num=-1):
        """Compute a Gaussian noise vector from the flux variance and the SNR, as computed from optical depth
        Parameters:
        snr : an array of signal to noise ratio (constant along each sightine)
        flux : an array of spectra (flux)  we want to add noise to
        spec_num : the index to spectra we want to add nose to. Leave it as -1 to add the noise to all spectra.
        """
        noise_array = np.array([])
        if np.size(np.shape(flux)) == 1:
            lines = 1
        else:
            lines = np.shape(flux)[0]
        #This is to get around the type rules.
        if lines == 1:
            #This ensures that we always get the same noise for the same spectrum
            np.random.seed(spec_num)
            flux += np.random.normal(0, 1./snr[spec_num], self.nbins)
        else:
            for ii in xrange(lines):
                np.random.seed(ii)
                noise = np.random.normal(0, 1./snr[ii], self.nbins)
                noise_array = np.append(noise_array, noise)
                flux[ii]+= noise
        return (flux, noise_array)


    def add_cont_error(self, CE, flux, spec_num=-1, u_delta=0.6, l_delta=-0.6):
        """Adding the Continuum error to spectra. If you want to add both random noise and continuum error, first add 
        the continuum error and then the random noise.
        Parameters:
        CE : the stdev of the gaussian noise 
        flux : an array of spectra (flux)  we want to add noise to
        spec_num : the index to spectra we want to add nose to. Leave it as -1 to add the noise to all spectra.
        u_delta, l_delta : upper and lower limit of the delta parameter
        """
        if np.size(np.shape(flux)) == 1:
            lines = 1
        else:
            lines = np.shape(flux)[0]
        #This is to get around the type rules
        if lines == 1:
            delta_array = np.array([])
            #This ensures that we always get the same noise for the same spectrum and is differen from seed for rand noise
            np.random.seed(2*spec_num)
            delta = np.random.normal(0, CE[spec_num])
            # Use lower and upper limit of delta from 2sigma for the highest CE in the survey
            while (delta < l_delta) or (delta > u_delta):
                delta = np.random.normal(0, CE[spec_num])
            flux /= (1.0 + delta)
        else:
            delta = np.empty(lines)
            for ii in xrange(lines):
                np.random.seed(2*ii)
                delta[ii] = np.random.normal(0, CE[ii])
                while (delta[ii] < l_delta) or (delta[ii] > u_delta):
                    delta[ii] = np.random.normal(0, CE[ii])
                flux[ii,:] /= (1.0 + delta[ii])        
        return (flux , delta)


    def load_savefile(self, savefile=None):
        """Load data from a file"""
        #Name of savefile
        try:
            f = h5py.File(savefile, 'r')
        except IOError:
            raise IOError("Could not read saved data from: "+savefile+". If the file does not exist, try using reload_file=True")
        grid_file = f["Header"]
        self.red = grid_file.attrs["redshift"]
        self.atime = 1./(1+self.red)
        self.OmegaM = grid_file.attrs["omegam"]
        self.nbins = grid_file.attrs["nbins"]
        self.omegab = grid_file.attrs["omegab"]
        self.OmegaLambda = grid_file.attrs["omegal"]
        self.hubble = grid_file.attrs["hubble"]
        self.npart = np.array(grid_file.attrs["npart"])
        self.box = grid_file.attrs["box"]
        self.discarded = grid_file.attrs["discarded"]
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
                    self.tau[(elem, int(ion), int(float(line)))] = np.array([0])
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
        try:
            grp = f["density_weight_density"]
            for elem in grp.keys():
                for ion in grp[elem].keys():
                    self.dens_weight_dens[(elem, int(ion))] = np.array([0])
        except KeyError:
            pass
        grp = f["num_important"]
        for elem in grp.keys():
            for ion in grp[elem].keys():
                self.num_important[(elem, int(ion))] = np.array(grp[elem][ion])
        grp = f["spectra"]
        self.cofm = np.array(grp["cofm"])
        self.axis = np.array(grp["axis"])
        f.close()

    def _interpolate_single_file(self, nsegment, elem, ion, ll, get_tau, load_all_data_first=False):
        """Read arrays and perform interpolation for a single file"""
        (pos, vel, elem_den, temp, hh, amumass) = self._read_particle_data(nsegment, elem, ion, get_tau)
        if load_all_data_first:
            for nseg in range(1, self.snapshot_set.get_n_segments()):
                (pos_, vel_, elem_den_, temp_, hh_, amumass_) = self._read_particle_data(nseg, elem, ion, get_tau)
                if amumass_ is False:
                    continue
                # We cannot concatenate onto empty arrays,
                #so if the first segment contained no particles we must rename
                if amumass is False:
                    pos = pos_
                    vel = vel_
                    temp = temp_
                    elem_den = elem_den_
                    hh = hh_
                else:
                    pos = np.concatenate((pos, pos_), axis=0)
                    if get_tau:
                        vel = np.concatenate((vel, vel_), axis=0)
                    elem_den = np.append(elem_den, elem_den_)
                    if get_tau or (ion != -1 and elem != 'H'):
                        temp = np.append(temp, temp_)
                    hh = np.append(hh, hh_)
                amumass = amumass_
        if amumass is False:
            return np.zeros([np.shape(self.cofm)[0], self.nbins], dtype=np.float32)
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
                line = self.lines[(elem, ion)][ll]
        else:
            #If get_tau is false, we don't use the line data
            #(only the density), so we can set this as we like.
            #Setting something makes writing the interfaces easier,
            #because C doesn't have default arguments.
            line = self.lines[("H", 1)][1215]
        # print(pos.shape, vel.shape, elem_den.shape, temp.shape, hh.shape, amumass)
        return self._do_interpolation_work(pos, vel, elem_den, temp, hh, amumass, line, get_tau)

    def _read_particle_data(self, fn, elem, ion, get_tau):
        """Read the particle data for a single interpolation"""
        pos = self.snapshot_set.get_data(0, "Position", segment=fn).astype(np.float32)
        hh = self.snapshot_set.get_smooth_length(0, segment=fn).astype(np.float32)

        #Find particles we care about
        if self.cofm_final:
            try:
                ind = self.part_ind[fn]
            except KeyError:
                ind = self.particles_near_lines(pos, hh, self.axis, self.cofm)
                self.part_ind[fn] = ind
        else:
            ind = self.particles_near_lines(pos, hh, self.axis, self.cofm)
        #Do nothing if there aren't any, and return a suitably shaped zero array
        if np.size(ind) == 0:
            return (False, False, False, False, False, False)
        pos = pos[ind, :]
        hh = hh[ind]
        #Get the rest of the arrays: reducing them each time to have a smaller memory footprint
        vel = np.zeros(1, dtype=np.float32)
        temp = np.zeros(1, dtype=np.float32)
        if get_tau:
            vel = self.snapshot_set.get_peculiar_velocity(0, segment=fn).astype(np.float32)
            vel = vel[ind, :]
        #gas density amu / cm^3
        den = self.gasprop.get_code_rhoH(0, segment=fn).astype(np.float32)
        # Get mass of atomic species
        if elem != "Z":
            amumass = self.lines.get_mass(elem)
        else:
            amumass = 1
        den = den[ind]
        #Only need temp for ionic density, and tau later
        if get_tau or (ion != -1 and elem != 'H'):
            temp = self.snapshot_set.get_temp(0, segment=fn, units=self.units).astype(np.float32)
            temp = temp[ind]
            #Some codes occasionally output negative temperatures, fix them
            it = np.where(temp <= 0)
            temp[it] = 1
        #Find the mass fraction in this ion
        #Get the mass fraction in this species: elem_den is now density in ionic species in amu/cm^2 kpc/h
        #(these weird units are chosen to be correct when multiplied by the smoothing length)
        elem_den = (den*self.rscale)*self.get_mass_frac(elem, fn, ind)
        #Special case H1:
        if elem == 'H' and ion == 1:
            # Neutral hydrogen mass frac
            elem_den *= (self.gasprop.get_reproc_HI(0, segment=fn)[ind]).astype(np.float32)
        elif ion != -1:
            #Cloudy density in physical H atoms / cm^3
            ind2 = self._filter_particles(elem_den, pos, vel, den)
            if np.size(ind2) == 0:
                return (False, False, False, False, False, False)
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
        elem_den /= amumass
        #Do interpolation.
        return (pos, vel, elem_den, temp, hh, amumass)

    def find_all_particles(self):
        """Returns the positions, velocities and smoothing lengths of all particles near sightlines."""
        nsegments = self.snapshot_set.get_n_segments()
        pp = np.empty([0, 3])
        hhh = np.array([])
        for i in range(nsegments):
            (pos, _, _, _, hh, amumass) = self._read_particle_data(i, "H", -1, False)
            if amumass is not False:
                pp = np.concatenate([pp, pos])
                hhh = np.concatenate([hhh, hh])
        return pp, hhh

    def _filter_particles(self, elem_den, pos, velocity, den):
        """Get a filtered list of particles to add to the sightlines"""
        _ = (pos, velocity, den)
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

    def _do_interpolation_work(self, pos, vel, elem_den, temp, hh, amumass, line, get_tau):
        """Run the interpolation on some pre-determined arrays, spat out by _read_particle_data"""
        #Factor of 10^-8 converts line width (lambda_X) from Angstrom to cm
        return _Particle_Interpolate(get_tau*1, self.nbins, self.kernel_int, self.box, self.velfac, self.atime, line.lambda_X*1e-8, line.gamma_X, line.fosc_X, amumass, self.tautail, pos, vel, elem_den, temp, hh, self.axis, self.cofm)

    def particles_near_lines(self, pos, hh, axis=None, cofm=None):
        """Filter a particle list, returning an index list of those near sightlines"""
        if axis is None:
            axis = self.axis
        if cofm is None:
            cofm = self.cofm
        #Axis is 1-indexed between 1 and 3. 1 is x axis.
        assert np.min(axis) > 0
        assert np.max(axis) < 4
        ind = _near_lines(self.box, pos, hh, axis, cofm)
        return ind

    def get_mass_frac(self, elem, fn, ind):
        """Get the mass fraction of a given species from a snapshot.
        Arguments:
            elem = name of element
            data = pointer to hdf5 array containing baryons
            ind = index of particles we care about
        Returns mass_frac - mass fraction of this ion
        """
        if elem == "Z":
            mass_frac = self.snapshot_set.get_data(0, "Metallicity", segment=fn).astype(np.float32)
        else:
            nelem = self.species.index(elem)
            #Get metallicity of this metal species
            try:
                mass_frac = (self.snapshot_set.get_data(0, "GFM_Metals", segment=fn).astype(np.float32))[:, nelem]
            except KeyError:
                #If GFM_Metals is not defined, fall back to primordial abundances
                metal_abund = np.array([0.76, 0.24], dtype=np.float32)
                nvalues = self.snapshot_set.get_blocklen(0, "Density", segment=fn)
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
        col_den = self.compute_spectra(elem, ion, 1215, False)
        ind = self.filter_DLA(col_den, thresh)
        H1_DLA = np.empty_like(col_den)
        #Update saves
        top = np.min([wanted, found+np.size(ind)])
        cofm_DLA[found:top] = self.cofm[ind][:top, :]
        H1_DLA[found:top] = col_den[ind][:top, :]
        found += np.size(ind)
        self.discarded = self.NumLos-np.size(ind)
        print("Discarded: ", self.discarded)
        while found < wanted:
            #Get a bunch of new spectra
            self.cofm = self.get_cofm()
            col_den = self.compute_spectra(elem, ion, 1215, False)
            ind = self.filter_DLA(col_den, thresh)
            #Update saves
            top = np.min([wanted, found+np.size(ind)])
            cofm_DLA[found:top] = self.cofm[ind][:top-found, :]
            H1_DLA[found:top] = col_den[ind][:top-found, :]
            found += np.size(ind)
            self.discarded += self.NumLos-np.size(ind)
            print("Discarded: ", self.discarded)
        #Correct proportions in case we find slightly more than we need
        self.discarded = int(self.discarded*1.*wanted/1./found)
        #Copy back
        self.cofm = cofm_DLA
        self.axis = self.axis[:ndla]
        self.colden[("H", 1)] = H1_DLA[:top]
        #Finalise the cofm array
        self.cofm_final = True
        self.NumLos = ndla

    def get_cofm(self, num=None):
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
        MM = self.get_density("Z", -1)
        HH = self.get_density("H", -1)
        if width > 0:
            (roll, hhr) = spec_utils.get_rolled_spectra(HH)
            mmr = np.array([np.roll(MMr, rr) for (MMr, rr) in zip(MM, roll)])
            imax = int(np.shape(MM)[1]/2)
            wbins = min(int(width/self.dvbin), imax)
            mms = np.array([np.sum(mmrr[imax-wbins:imax+wbins]) for mmrr in mmr])
            hhs = np.array([np.sum(hhrr[imax-wbins:imax+wbins]) for hhrr in hhr])
        else:
            mms = np.sum(MM, axis=1)
            hhs = np.sum(HH, axis=1)
        return mms/hhs/self.solarz
        #Use only DLA regions: tricky to preserve shape
        #ma_HH = np.ma.masked_where(HH < thresh, MM/HH)
        #data = np.array([np.mean(ma_HH, axis=1)])
        #return data/solar

    def get_ion_metallicity(self, species, ion):
        """Get the metallicity derived from an ionic species"""
        MM = self.get_density(species, ion)
        HH = self.get_density("H", 1)
        mms = np.sum(MM, axis=1)
        hhs = np.sum(HH, axis=1)
        return mms/hhs/self.solar[species]

    def compute_spectra(self, elem, ion, ll, get_tau):
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
        nsegments = self.snapshot_set.get_n_segments(part_type=0)
        arepo = (self.kernel_int == 2)
        if arepo :
           nsegments=1
        result = self._interpolate_single_file(0, elem, ion, ll, get_tau, load_all_data_first=arepo)
        #Do remaining files
        for nn in xrange(1, nsegments):
            tresult = self._interpolate_single_file(nn, elem, ion, ll, get_tau)
            print("Interpolation %.1f percent done" % (100*nn/nsegments), flush=True)
            #Add new file
            result += tresult
            del tresult
        if self.MPI is None :
           return result
        else :
           # Make sure the data is contiguous in memory
           result = np.ascontiguousarray(result, np.float32)
           # Each rank constructs a portion of the spectrum. Add all the portions
           self.comm.Allreduce(self.MPI.IN_PLACE, result, op=self.MPI.SUM)
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
        eq_width = np.trapz(-np.expm1(-tau), dx=dl, axis=1)
        #Don't need to divide by 1+z as lambda_X is already rest wavelength
        return eq_width

    def eq_width_hist(self, elem, ion, line, dv=0.05):
        """
        Compute a histogram of the equivalent width distribution of our spectra.

        Returns:
            (v, f_table) - v (binned in log) and corresponding f(N)
        """
        print("For ", line, " Angstrom")
        eq_width = self.equivalent_width(elem, ion, line)
        ii = np.where(eq_width > 0)
        v_table = np.arange(np.log10(np.min(eq_width[ii])), np.log10(np.max(eq_width[ii])), dv)
        vbin = (v_table[1:]+v_table[:-1])/2.
        eqws = np.histogram(np.log10(eq_width[ii]), v_table, density=True)[0]
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
            colden = self.compute_spectra(elem, ion, 0, False)
            self.colden[(elem, ion)] = colden
            return colden

    def get_density(self, elem, ion, force_recompute=False):
        """Get the density in each pixel for a given species, in units of [metal] ions cm^-3."""
        colden = self.get_col_density(elem, ion, force_recompute)
        phys = self.dvbin/self.velfac*self.rscale
        return colden/phys

    def get_tau(self, elem, ion, line, number=-1, force_recompute=False):
        """Get the optical depth in each pixel along the sightline for a given line."""
        try:
            if force_recompute:
                raise KeyError
            self._really_load_array((elem, ion, line), self.tau, "tau")
            tau = self.tau[(elem, ion, line)]
        except KeyError:
            tau = self.compute_spectra(elem, ion, line, True)
            self.tau[(elem, ion, line)] = tau
        if number >= 0:
            tau = tau[number, :]
        return tau


    def get_observer_tau(self, elem, ion, number=-1, force_recompute=False, noise=True):
        """Get the optical depth for a particular element out of:
           (He, C, N, O, Ne, Mg, Si, Fe)
           and some ion number, choosing the line which causes the maximum optical depth to be closest to unity.
        """
        try:
            if force_recompute:
                raise KeyError
            self._really_load_array((elem, ion), self.tau_obs, "tau_obs")
            ntau = self.tau_obs[(elem, ion)]
        except KeyError:
            #Compute tau for each line
            nlines = len(self.lines[(elem, ion)])
            tau = np.zeros([nlines, self.NumLos, self.nbins])
            for ll in range(nlines):
                line = list(self.lines[(elem, ion)].keys())[ll]
                tau_loc = self.compute_spectra(elem, ion, line, True)
                tau[ll, :, :] = tau_loc
                del tau_loc
            #Maximum tau in each spectra with each line,
            #after convolving with a Gaussian for instrumental broadening.
            maxtaus = np.max(spec_utils.res_corr(tau, self.dvbin, self.spec_res), axis=-1)
            #Array for line indices
            ntau = np.empty([self.NumLos, self.nbins])
            #Use the maximum unsaturated optical depth
            for ii in xrange(self.NumLos):
                # we want unsaturated lines, defined as those with tau < 3
                #which is the maximum tau in the sample of Neeleman 2013
                #Also use lines with some absorption: tau > 0.1, roughly twice noise level.
                ind = np.where(np.logical_and(maxtaus[:, ii] < 3, maxtaus[:, ii] > 0.1))
                if np.size(ind) > 0:
                    line = np.where(maxtaus[:, ii] == np.max(maxtaus[ind, ii]))
                else:
                    #We have no lines in the desired region: here use something slightly saturated.
                    #In reality the observers will use a different ion
                    ind2 = np.where(maxtaus[:, ii] > 0.1)
                    if np.size(ind2) > 0:
                        line = np.where(maxtaus[:, ii] == np.min(maxtaus[ind2, ii]))
                    else:
                        #We have no observable lines: this spectra are metal-poor
                        #and will be filtered anyway.
                        line = np.where(maxtaus[:, ii] == np.max(maxtaus[:, ii]))
                if np.size(line) > 1:
                    line = (line[0][0], )
                ntau[ii, :] = tau[line, ii, :]
            self.tau_obs[(elem, ion)] = ntau
        if number >= 0:
            ntau = ntau[number, :]
        # Convolve lines by a Gaussian filter of the resolution of the spectrograph.
        ntau = spec_utils.res_corr(ntau, self.dvbin, self.spec_res)
        #Add noise
        if noise and self.snr > 0:
            ntau = self.add_noise(self.snr, ntau, number)
        return ntau

    def vel_width(self, elem, ion):
        """
           Find the velocity width of an ion.
           This is the width in velocity space containing 90% of the optical depth
           over the absorber.
           elem - element to look at
           ion - ionisation state of this element.
        """
        try:
            return self.vel_widths[(elem, ion)]
        except KeyError:
            tau = self.get_observer_tau(elem, ion)
            (low, high, offset) = self.find_absorber_width(elem, ion)
            #  Size of a single velocity bin
            vel_width = np.zeros(np.shape(tau)[0])
            #deal with periodicity by making sure the deepest point is in the middle
            for ll in np.arange(0, np.shape(tau)[0]):
                tau_l = np.roll(tau[ll, :], offset[ll])[low[ll]:high[ll]]
                (nnlow, nnhigh) = self._vel_width_bound(tau_l)
                vel_width[ll] = self.dvbin*(nnhigh-nnlow)
            #Return the width
            self.vel_widths[(elem, ion)] = vel_width
            return self.vel_widths[(elem, ion)]


    def _vel_single_file(self, fn, elem, ion):
        """Get the column density weighted interpolated velocity field for a single file"""
        (pos, vel, elem_den, temp, hh, amumass) = self._read_particle_data(fn, elem, ion, True)
        if amumass is False:
            return np.zeros([np.shape(self.cofm)[0], self.nbins, 3], dtype=np.float32)
        line = self.lines[("H", 1)][1215]
        vv = np.empty([np.shape(self.cofm)[0], self.nbins, 3], dtype=np.float32)
        phys = self.dvbin/self.velfac*self.rscale
        for ax in (0, 1, 2):
            weight = vel[:, ax]*np.sqrt(self.atime)
            vv[:, :, ax] = self._do_interpolation_work(pos, vel, elem_den*weight/phys, temp, hh, amumass, line, False)
        return vv

    def _get_mass_weight_quantity(self, func, elem, ion):
        """
        Helper function to get a mass weighted quantity, which reduces code duplication.
        func should be something like _vel_single_file above (for velocity)
        and have the signature func(file, elem, ion)
        """
        nsegments = self.snapshot_set.get_n_segments()
        result = func(0, elem, ion)
        if nsegments > 1:
            #Do remaining files
            for nn in xrange(1, nsegments):
                tresult = func(nn, elem, ion)
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
                result[:, :, ax] /= den
        return result

    def get_velocity(self, elem, ion):
        """Get the column density weighted velocity in each pixel for a given species.
        """
        try:
            self._really_load_array((elem, ion), self.velocity, "velocity")
            return self.velocity[(elem, ion)]
        except KeyError:
            velocity = self._get_mass_weight_quantity(self._vel_single_file, elem, ion)
            self.velocity[(elem, ion)] = velocity
            return velocity

    def _temp_single_file(self, fn, elem, ion):
        """Get the density weighted interpolated temperature field for a single file"""
        (pos, vel, elem_den, temp, hh, amumass) = self._read_particle_data(fn, elem, ion, True)
        if amumass is False:
            return np.zeros([np.shape(self.cofm)[0], self.nbins], dtype=np.float32)
        line = self.lines[("H", 1)][1215]
        phys = np.float32(self.dvbin/self.velfac*self.rscale)
        temp = self._do_interpolation_work(pos, vel, elem_den*temp/phys, temp, hh, amumass, line, False)
        return temp

    def get_temp(self, elem, ion):
        """Get the density weighted temperature in each pixel for a given species.
        """
        try:
            self._really_load_array((elem, ion), self.temp, "temperature")
            return self.temp[(elem, ion)]
        except KeyError:
            temp = self._get_mass_weight_quantity(self._temp_single_file, elem, ion)
            self.temp[(elem, ion)] = temp
            return temp

    def _densweightdens(self, fn, elem, ion):
        """Get the density weighted interpolated density field for a single file"""
        (pos, vel, elem_den, temp, hh, amumass) = self._read_particle_data(fn, elem, ion, True)
        if amumass is False:
            return np.zeros([np.shape(self.cofm)[0], self.nbins], dtype=np.float32)
        (_, _, species_den, _, _, _) = self._read_particle_data(fn, elem, -1, True)
        line = self.lines[("H", 1)][1215]
        phys = np.float32(self.dvbin/self.velfac*self.rscale)
        dens = self._do_interpolation_work(pos, vel, (elem_den/phys)*(species_den/self.rscale), temp, hh, amumass, line, False)
        return dens

    def get_dens_weighted_density(self, elem, ion):
        """This function gets the (ion) density weighted (species) density in each pixel for a given species.
        This may seem an odd thing to compute, but it shows the characteristic density which dominates the
        production of a given ionic species is found in these spectra.

        For example, to see which gas densities produce the Lyman alpha forest, one could do:
            get_dens_weighted_density("H", 1)
        which would be:
            int ( rho_HI rho dv) / int(rho_HI dv)
        where rho is the density of hydrogen and rho_HI is the density of neutral hydrogen.
        """
        try:
            self._really_load_array((elem, ion), self.dens_weight_dens, "density_weight_density")
            return self.dens_weight_dens[(elem, ion)]
        except KeyError:
            dens_weight_dens = self._get_mass_weight_quantity(self._densweightdens, elem, ion)
            self.dens_weight_dens[(elem, ion)] = dens_weight_dens
            return dens_weight_dens

    def get_b_param_dist(self, elem="H", ion=1, line=1215):
        """Get the power law betweeen the 'minimum' b parameter and column density."""
        tau = self.get_tau(elem, ion, line)
        (b_0, gamm1) = voigtfit.get_b_param_dist(tau, self.dvbin, elem=elem, ion=ion, line=line)
        return (b_0, gamm1)

    def column_density_from_voigt(self, elem="H", ion=1, line=1215, dlogN=0.2, minN=13, maxN=23, close=50., dX=True, nspectra=-1):
        """This computes the column density function using column densities from a Voigt profile fit.
        Concatenate objects closer than close km/s."""
        NHI_table = 10**np.arange(minN, maxN, dlogN)
        center = np.array([(NHI_table[i]+NHI_table[i+1])/2. for i in range(0, np.size(NHI_table)-1)])
        width = np.array([NHI_table[i+1]-NHI_table[i] for i in range(0, np.size(NHI_table)-1)])
        #Number of lines
        tot_lines = self.NumLos+self.discarded
        #Absorption distance for each line
        if dX:
            dist = self.units.absorption_distance(self.box, self.red)
        else:
            dist = self.units.redshift_distance(self.box, self.red, self.OmegaM)
        tau = self.get_tau(elem, ion, line)
        if nspectra > 0:
            tau = tau[:nspectra, :]
            tot_lines *= nspectra/self.NumLos
        #Combine all lines closer than the close value into one.
        n_vals = voigtfit.get_voigt_systems(tau, self.dvbin, elem=elem, ion=ion, line=line, verbose=False, close=close)
        (tot_f_N, NHI_table) = np.histogram(n_vals, NHI_table)
        #The normalisation should be the total sightline distance.
        tot_f_N = tot_f_N/(width*dist*tot_lines)
        return (center, tot_f_N)

    def column_density_function(self, elem="H", ion=1, dlogN=0.2, minN=13, maxN=23., line=True, close=50., dX=True):
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
        center = np.array([(NHI_table[i]+NHI_table[i+1])/2. for i in range(0, np.size(NHI_table)-1)])
        width = np.array([NHI_table[i+1]-NHI_table[i] for i in range(0, np.size(NHI_table)-1)])
        #Number of lines
        tot_lines = self.NumLos+self.discarded
        #Absorption distance for each line
        if dX:
            dist = self.units.absorption_distance(self.box, self.red)
        else:
            dist = self.units.redshift_distance(self.box, self.red, self.OmegaM)
        if line:
            rho = np.sum(self.get_col_density(elem, ion), axis=1)
        else:
            rho = self.get_col_density(elem, ion)
            cbins = np.max((int(np.round((close/self.dvbin))), 1))
            rhob = np.array([np.sum(rho[:, cbins*i:cbins*(i+1)], axis=1) for i in xrange(int(np.shape(rho)[1]/cbins))]).T
            #Check that fp roundoff is not too severe: this can sometimes trigger for silly reasons
            #assert np.abs((np.sum(rhob) / np.sum(rho))-1) < 5e-2
            rho = rhob
        (tot_f_N, NHI_table) = np.histogram(rho, NHI_table)
        #The normalisation should be the total sightline distance.
        tot_f_N = tot_f_N/(width*dist*tot_lines)
        return (center, tot_f_N)

    def _rho_abs(self, thresh=10**20.3, upthresh=None, elem="H", ion=1):
        """Compute rho_abs, the sum of the mass in an absorber,
           divided by the volume of the spectra in g/cm^3 (comoving).
            Omega_DLA = m_p * avg. column density / (1+z)^2 / length of column
        """
        #Column density of ion in atoms cm^-2 (physical)
        col_den = np.sum(self.get_col_density(elem, ion), axis=1)
        if thresh > 0 or upthresh is not None:
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

    def omega_abs(self, thresh=10**20.3, upthresh=1e40, elem="H", ion=1):
        """Compute Omega_abs, the sum of the mass in a given absorber,
            divided by the volume of the spectra, divided by the critical density.
            Ω_abs = m_p * avg. column density / (1+z)^2 / length of column / rho_c
            Note: If we want the neutral gas density rather than the neutral hydrogen density, divide by 0.76,
            the hydrogen mass fraction.
        """
        #Avg density in g/cm^3 (comoving) divided by critical density in g/cm^3
        omega_DLA = self._rho_abs(thresh, upthresh, elem=elem, ion=ion)/self.units.rho_crit(self.hubble)
        return omega_DLA

    def omega_abs_cddf(self, thresh=10**20.3, upthresh=1e40, elem="H", ion=1):
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
        h100 = self.units.h100*self.hubble
        #The 1+z factor converts lightspeed to comoving
        omega_abs = self.lines.get_mass(elem)*self.units.protonmass/self.units.light*h100/self.units.rho_crit(self.hubble)*np.trapz(moment, bins)
        return omega_abs

    def line_density(self, thresh=10**20.3, upthresh=10**40, elem="H", ion=1):
        """Compute the line density, the total no. of DLA sightlines divided by the total number of sightlines, multiplied by d L / dX. This is dN/dX = l_DLA(z)
        """
        col_den = np.sum(self.get_col_density(elem, ion), axis=1)
        #Average fraction of pixels containing a DLA
        frac = 1.*np.size(col_den[np.where((col_den > thresh)*(col_den < upthresh))])/np.size(col_den)
        #Divide by abs. distance per sightline
        frac *= np.size(col_den)/(np.size(col_den)+1.*self.discarded)
        return frac/(self.units.absorption_distance(self.box, self.red))

    def line_density_eq_w(self, thresh=0.4, elem="H", ion=1, line=1216):
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
            spos = cofm[:, 1:]
        if axis == 2:
            spos = np.vstack([cofm[:, 0], cofm[:, 2]]).T
        if axis == 3:
            spos = cofm[:, :2]
        return spos

    def _filter_tau(self, tau, tau_thresh=None):
        """Filter optical depths to remove sightlines with optically thick absorbers."""
        if tau_thresh is not None:
            tausum = np.max(tau, axis=1)
            ii = np.where(tausum < tau_thresh)
            tau = tau[ii]
        return tau

    def get_mean_flux(self, elem="H", ion=1, line=1215, tau_thresh=None):
        """Get the mean flux along a set of sightlines"""
        tau = self.get_tau(elem, ion, line)
        tau = self._filter_tau(tau, tau_thresh=tau_thresh)
        return np.mean(np.exp(-tau))

    def get_flux_pdf(self, elem="H", ion=1, line=1215, nbins=20, mean_flux_desired=None, tau_thresh=None):
        """Get the flux PDF, a histogram of the flux values."""
        tau = self.get_tau(elem, ion, line)
        tau = self._filter_tau(tau, tau_thresh=tau_thresh)
        return fstat.flux_pdf(tau, nbins=nbins, mean_flux_desired=mean_flux_desired)

    def get_flux_power_1D(self, elem="H", ion=1, line=1215, mean_flux_desired=None, window=True, tau_thresh=None):
        """Get the power spectrum of (variations in) the flux along the line of sight.
        This is: P_F(k_F) = <d_F d_F>
                 d_F = e^-tau / mean(e^-tau) - 1
        Arguments:
            mean_flux_desired: if not None, the spectral optical depths will be rescaled
                to match the desired mean flux.
            window: if True, the flux power spectrum is divided by the window function for the pixel width.
                    This interacts poorly with mean flux rescaling.
            tau_thresh: sightlines with a total optical depth greater than this value are removed before mean flux rescaling."""
        tau = self.get_tau(elem, ion, line)
        #Remove sightlines which contain a strong absorber
        tau = self._filter_tau(tau, tau_thresh=tau_thresh)
        #Mean flux rescaling does not commute with the spectrum resolution correction!
        if mean_flux_desired is not None and window is True and self.spec_res > 0:
            raise ValueError("Cannot sensibly rescale mean flux with gaussian smoothing")
        (kf, avg_flux_power) = fstat.flux_power(tau, self.vmax, spec_res=self.spec_res, mean_flux_desired=mean_flux_desired, window=window)
        return kf[1:], avg_flux_power[1:]

    def spline_fit(self, flux_i, chi_min=3., vel_seg_min=10., ini_break_spacing=50.):
        """
        Fit flux spectra (converge chi^2) using a cubic spline with adaptive breakpoints.

        chi_min = minimum change in chi^2 to accept a new breakpoint
        vel_seg_min = segment size limit (km/sec) i.e. resolution minimum
        ini_break_spacing = initial breakpoint spacing (km/sec)

        Returns: splines for each sightline (tau.shape)
        Note -- splines may return flux outside the 0 < F < 1 range
        """
        from scipy.interpolate import LSQUnivariateSpline
        def chi_squared(expected, observed):
        # observed is the fit, expected is the original spectrum
            return np.sum(((expected - observed)/np.std(expected))**2)

        # array of velocities (km/sec), i.e. x-axis
        vel = np.linspace(0, self.vmax, flux_i.shape[1])
        vel_stepsize = vel[1]-vel[0] # velocity bin size (km/sec)
        if vel_stepsize >= vel_seg_min:
            raise Exception("Velocity resolution must be less than minimum segment size (vel_stepsize < vel_seg_min)")

        # index spacing to get ~velocity breakpoint spacing
        ind_break_spacing = int(np.round(ini_break_spacing/vel_stepsize))

        all_spline_flux = np.zeros(flux_i.shape)
        for j in range(self.NumLos):

            flux = flux_i[j] # flux for current sight line

            # initial interior breakpoints
            interior_inds = np.arange(ind_break_spacing, vel.size, ind_break_spacing)
            # ensure the last segment is larger than vel_seg_min
            # 'segments' are the data between pairs of breakpoints
            if (vel[-1] - vel[interior_inds[-1]]) < vel_seg_min:
                interior_inds = interior_inds[:-1]

            # initialize number of unconverged segments (to start while loop)
            n_unconverged = 1
            while n_unconverged > 0: # while there are unconverged segments
                # fit the spline using the current step breakpoints
                current_spline = LSQUnivariateSpline(vel, flux, vel[interior_inds], k=3, ext=3)

                # indices for breakpoints, with start and end (interior + boundary points)
                all_inds = np.concatenate([[0], interior_inds, [vel.size]])
                # indices for breakpoints for maximal next step
                # (i.e. with breakpoints added between current pairs)
                next_inds = np.sort(np.concatenate([all_inds, (all_inds[:-1] + all_inds[1:])//2]))[1:-1]
                # fit the spline using the next step breakpoints, next_inds
                next_spline = LSQUnivariateSpline(vel, flux, vel[next_inds], k=3, ext=3)

                # initialize array to contain new breakpoints to be added
                new_bp = np.empty(0, dtype=int)
                for i in range(all_inds.size-1): # for each segment
                    # current segment velocities
                    seg_vel = vel[all_inds[i]:all_inds[i+1]+1]
                    # current segment flux
                    seg_flux = flux[all_inds[i]:all_inds[i+1]+1]
                    # current segment spline flux
                    seg_spline_flux = current_spline(seg_vel)

                    # chi squared for current spline versus input flux
                    chisq = chi_squared(seg_flux, seg_spline_flux)
                    # chisq for SAME segments, with added breakpoints
                    next_chisq = chi_squared(seg_flux, next_spline(seg_vel))

                    # segments sizes that would result from including proposed breakpoint
                    new_segs = np.array([seg_vel.max()-vel[next_inds[2*i]], vel[next_inds[2*i]]-seg_vel.min()])
                    # if chisq changes more than chi_min AND new segments > vel_seg_min
                    if (chisq-next_chisq) > chi_min and new_segs[0] >= vel_seg_min and new_segs[1] >= vel_seg_min:
                        # add a breakpoint ~halfway in the segment (nearest index)
                        new_bp = np.append(new_bp, int(np.round(all_inds[i]/2 + all_inds[i+1]/2)))

                # number of segments that were improved with additional breakpoints
                n_unconverged = new_bp.size

                if n_unconverged > 0: # if there is at least 1 breakpoint to add
                    # add that breakpoint to the list of breakpoints and sort
                    interior_inds = np.sort(np.append(interior_inds, new_bp))

            # add current sight line spline to array of splines
            all_spline_flux[j] = current_spline(vel)

        return all_spline_flux


    def renormalize_flux(self, flux, section_size):
        """
        Renormalize flux into sections approximately section_size.

        flux = input flux, np.exp(-tau)
        section_size = desired section size

        Return: array of renormalized flux, arranged by section
        (i.e. the first section for all sight lines, then the 2nd section . . .)
        """

        # velocity step size (resolution, km/sec)
        vel_stepsize = self.vmax/(flux.shape[1]-1)

        if section_size < vel_stepsize/self.velfac or section_size > self.box:
            raise Exception("Section size must be greater than spatial resolution and <= box size.\n"+"Spatial resolution is "+str(vel_stepsize/self.velfac)+", box size is "+str(self.box))

        # number of (non-integer) sections that could fit into the box
        n_sections = self.box/section_size
        # number of indices for each section to get ~section_size
        n_inds = int(np.floor(self.vmax/(n_sections*vel_stepsize)))

        # initialize array for storing section flux
        # shape is [num_sight_lines * n_sections_per_sight_line, section length + 1]
        flux_sections = np.zeros([flux.shape[0]*int(np.floor(n_sections)), n_inds+1])

        # running counter for the flux_sections array index (where to store)
        fs_ind = 0
        for i in range(int(np.floor(n_sections))): # for each section

            for j in range(flux.shape[0]): # for each line of sight

                # store the ith section flux into the flux_sections array
                flux_sections[fs_ind] = flux[j, i*n_inds:(i+1)*n_inds+1]
                # and divide by the maximum in the section
                flux_sections[fs_ind] = flux_sections[fs_ind]/flux_sections[fs_ind].max()
                fs_ind = fs_ind + 1

        return flux_sections


    def compute_curvature(self, flux):
        """
        Compute mean absolute curvatures for a set of flux spectra.

        Returns: the mean absolute curvature for each spectra passed.
        """

        # velocities, needed for more accurate derivatives
        vel = np.linspace(0, self.vmax, self.nbins)[0:flux.shape[1]]
        # initialize an array for the mean absolute curvatures
        curvature = np.zeros(flux.shape[0])

        for i in range(curvature.size): # for each section

            # limit flux to 0.1 < F < 0.9
            flux_range = np.where(np.logical_and(flux[i] > 0.1, flux[i] < 0.9))[0]
            # check that enough values actually fall in this range for a second derivative
            if flux_range.size > 2:
                # get first derivative
                flux_p = np.gradient(flux[i, flux_range], vel[flux_range], edge_order=2)
                # and get second derivative
                flux_pp = np.gradient(flux_p, vel[flux_range], edge_order=2)

                # compute the mean absolute curvature
                curvature[i] = np.mean(np.abs(flux_pp/((1 + flux_p**2)**(3/2))))

        # ensure only sections with non-zero curvature are returned
        curvature = curvature[np.where(curvature > 0)]

        return curvature


    def get_curvature(self, elem, ion, line, section_size=0, snr_input=None, chi=3., seg_res=10., ini_break=50.):
        """
        Calculate the mean absolute curvature for a set of sight lines. Renormalizes each sight line spectra
        into sections of size section_size, then rescales them before computing the curvature.

        section_size = desired (spatial) size of spectra sections (spatial res <= section_size <= box size)
        snr_input = snr for the flux (array of size n_sightlines)
        chi = minimum chi^2 change to accept a new breakpoint in spline fitting
        seg_res = segment size limit (km/sec) i.e. resolution minimum for spline fitting
        ini_break = initial breakpoint spacing (km/sec) for spline fitting

        Returns: mean absolute curvature for each section of each sight line.
        """

        # check for noise in spectra, get spline fits accordingly
        # if no noise is requested, just return the flux
        if snr_input is None:
            flux = np.exp(-self.get_tau(elem, ion, line))

        # if a noise input is specified here, add that, then fit spline
        else: 
            flux_i = self.add_noise(snr_input, np.exp(-self.get_tau(elem, ion, line)))[0]
            flux = self.spline_fit(flux_i, chi_min=chi, vel_seg_min=seg_res, ini_break_spacing=ini_break)

        if section_size == 0: # if no section size is passed, assume boxsize
            section_size = self.box
        # renormalize in section_size chunks (i.e. divide spectra and normalize)
        flux_sections = self.renormalize_flux(flux, section_size)

        # rescale the renormalized sections
        # tau should be nearly all the optical depths (minus where flux <= 0)
        tau = -np.log(flux_sections[np.where(flux_sections > 0)])
        # get the scaling factor
        scale = fstat.mean_flux(tau, np.exp(-fstat.obs_mean_tau(1/self.atime - 1)))
        # scale the positive flux
        # non-positive flux will be removed in the call to compute_curvature
        flux_sections[np.where(flux_sections > 0)] = flux_sections[np.where(flux_sections > 0)]**scale

        # compute mean absolute curvature for each rescaled, renormalized section
        curvature = self.compute_curvature(flux_sections)

        return curvature

