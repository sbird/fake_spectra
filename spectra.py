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

import numpy as np
import hsml
import math
import convert_cloudy
import cold_gas
import line_data
import h5py
import hdfsim
from scipy.interpolate import UnivariateSpline
from scipy.integrate import cumtrapz
from scipy.ndimage.filters import gaussian_filter1d
import os.path as path
import shutil
import subfindhdf
import numexpr as ne
from _spectra_priv import _Particle_Interpolate, _near_lines

class Spectra:
    """Class to interpolate particle densities along a line of sight and calculate their absorption
        Arguments:
            num - Snapshot number
            base - Name of base directory for snapshot
            cofm - table of los positions, as [n, 3] shape array.
            axis - axis along which to put the sightline
            res (optional) - Spectra pixel resolution in km/s
            snr - If nonzero, add noise for the requested signal to noise for the spectra, when loading from disc.
            spec_res - Resolution of the spectrograph. The spectra will be convolved with a Gaussian of this FWHM on loading from disc.
    """
    def __init__(self,num, base,cofm, axis, res=1., cdir=None, savefile="spectra.hdf5", savedir=None, reload_file=False, snr = 0., spec_res = 8):
        #Various physical constants
        #Internal gadget mass unit: 1e10 M_sun/h in g/h
        self.UnitMass_in_g=1.989e43
        #Internal gadget length unit: 1 kpc/h in cm/h
        self.UnitLength_in_cm=3.085678e21
        #Speed of light in cm/s
        self.light = 2.99e10
        #proton mass in g
        self.protonmass=1.67262178e-24
        #Newton's constant in cm^3/g/s^2
        self.gravcgs = 6.674e-8
        self.num = num
        self.base = base
        #Empty dictionary to add results to
        self.tau_obs = {}
        self.tau = {}
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
        self.npart = 512**3
        #If greater than zero, will add noise to spectra when they are loaded.
        self.snr = snr
        self.spec_res = spec_res
        #Minimum length of spectra within which to look at metal absorption (in km/s)
        self.minwidth = 500.
        try:
            self.files = hdfsim.get_all_files(num, base)
            self.files.reverse()
        except IOError:
            pass
        if savedir == None:
            savedir = path.join(base,"snapdir_"+str(num).rjust(3,'0'))
        self.savefile = path.join(savedir,savefile)
        #Snapshot data
        if reload_file:
            print "Reloading from snapshot (savefile: ",self.savefile," )"
            self.cofm = cofm
            self.axis = np.array(axis, dtype = np.int32)
            ff = h5py.File(self.files[0], "r")
            self.box = ff["Header"].attrs["BoxSize"]
            self.red = ff["Header"].attrs["Redshift"]
            self.atime = ff["Header"].attrs["Time"]
            self.hubble = ff["Header"].attrs["HubbleParam"]
            self.OmegaM = ff["Header"].attrs["Omega0"]
            self.OmegaLambda = ff["Header"].attrs["OmegaLambda"]
            self.npart=ff["Header"].attrs["NumPart_Total"]+2**32*ff["Header"].attrs["NumPart_Total_HighWord"]
            #Calculate omega_baryon (approximately)
            mass_dm = ff["Header"].attrs["MassTable"][1]*ff["Header"].attrs["NumPart_ThisFile"][1]
            mass_bar = np.sum(ff["PartType0"]["Masses"])
            self.omegab = mass_bar/(mass_bar+mass_dm)*self.OmegaM
            ff.close()
        else:
            self.load_savefile(self.savefile)
        # Conversion factors from internal units
        self.rscale = (self.UnitLength_in_cm*self.atime)/self.hubble    # convert length to physical cm
        #  Calculate the length scales to be used in the box: Hz in km/s/Mpc
        Hz = 100.0*self.hubble * np.sqrt(self.OmegaM/self.atime**3 + self.OmegaLambda)
        #Convert comoving kpc/h to physical km/s
        self.velfac = self.atime / self.hubble * Hz/1000
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
        #Generate cloudy tables
        if cdir != None:
            self.cloudy_table = convert_cloudy.CloudyTable(self.red, cdir)
        else:
            self.cloudy_table = convert_cloudy.CloudyTable(self.red)
        #Line data
        self.lines = line_data.LineData()
        print self.NumLos, " sightlines. resolution: ", self.dvbin, " z=", self.red
        #Try to load a halo catalogue
        self.load_halo()

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
            self._load_all_multihash(self.colden, "velocity")
        except IOError:
            pass
        try:
            if path.exists(self.savefile):
                shutil.move(self.savefile,self.savefile+".backup")
            f=h5py.File(self.savefile,'w')
        except IOError:
            try:
                f=h5py.File(self.savefile,'w')
            except IOError:
                raise IOError("Could not open ",self.savefile," for writing")
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
        for key in array.keys():
            self._really_load_array(key, array, array_name)

    def _save_multihash(self,save_array, grp):
        """Save an array using a tuple key, like save_array[(elem, ion, line)]
        to a hierarchy of hdf groups below grp"""
        for (key, value) in save_array.iteritems():
            #Create directory hierarchy recursively
            gg = grp
            for ii in xrange(np.size(key)-1):
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
        varnce = np.var(1-flux, axis=-1)
        #This is to get around the type rules.
        if np.size(varnce) == 1:
            #This ensures that we always get the same noise for the same spectrum
            np.random.seed(seed)
            flux += np.random.normal(0, np.sqrt(varnce/snr), self.nbins)
        else:
#             flux += np.array([np.random.normal(0, np.sqrt(vv/snr), self.nbins) if vv > 0 else np.zeros(self.nbins) for vv in varnce])
            for ii in xrange(np.size(varnce)):
                np.random.seed(ii)
                if varnce[ii] > 0:
                    flux[ii]+=np.random.normal(0,np.sqrt(varnce[ii]/snr), self.nbins)
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
            raise IOError("Could not open file: "+savefile)
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
                    self.tau[(elem, int(ion),int(line))] = np.array([0])
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

    def _interpolate_single_file(self,fn, elem, ion, ll, get_tau):
        """Read arrays and perform interpolation for a single file"""
        (pos, vel, elem_den, temp, hh, amumass) = self._read_particle_data(fn, elem, ion,get_tau)
        if amumass == False:
            return np.zeros([np.shape(self.cofm)[0],self.nbins],dtype=np.float32)
        if get_tau:
            line = self.lines[(elem,ion)][ll]
        else:
            line = self.lines[("H",1)][1215]
        return self._do_interpolation_work(pos, vel, elem_den, temp, hh, amumass, line, get_tau)

    def _read_particle_data(self,fn, elem, ion, get_tau):
        """Read the particle data for a single interpolation"""
        ff = h5py.File(fn, "r")
        data = ff["PartType0"]
        pos = np.array(data["Coordinates"],dtype=np.float32)
        hh = hsml.get_smooth_length(data)

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
        star=cold_gas.RahmatiRT(self.red, self.hubble)
        vel = np.zeros(1,dtype=np.float32)
        temp = np.zeros(1,dtype=np.float32)
        if get_tau:
            vel = np.array(data["Velocities"],dtype=np.float32)
            vel = vel[ind,:]
        #gas density amu / cm^3
        den=star.get_code_rhoH(data)
        # Get mass of atomic species
        if elem != "Z":
            amumass = self.lines.get_mass(elem)
        else:
            amumass = 1
        den = den[ind]
        #Only need temp for ionic density, and tau later
        if get_tau or (ion != -1 and elem != 'H'):
            temp = star.get_temp(data)
            temp = temp[ind]
        #Find the mass fraction in this ion
        #Get the mass fraction in this species: elem_den is now density in ionic species in amu/cm^2 kpc/h
        #(these weird units are chosen to be correct when multiplied by the smoothing length)
        elem_den = (den*self.rscale)*self.get_mass_frac(elem,data,ind)
        #Special case H1:
        if elem == 'H' and ion == 1:
            # Neutral hydrogen mass frac
            elem_den *= star.get_reproc_HI(data)[ind]
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
            elem_den = elem_den[ind2] * self._get_elem_den(elem, ion, den[ind2], temp, data, ind, ind2, star)
            del ind2
        ff.close()
        #Get rid of ind so we have some memory for the interpolator
        del den
        #Put density into number density of particles, from amu
        elem_den/=amumass
        #Do interpolation.
        return (pos, vel, elem_den, temp, hh, amumass)

    def _filter_particles(self, elem_den, pos, velocity, den):
        """Get a filtered list of particles to add to the sightlines"""
        ind2 = np.where(elem_den > 0)
        return ind2

    def _get_elem_den(self, elem, ion, den, temp, data, ind, ind2, star):
        """Get the density in an elemental species. Broken out so it can be over-ridden by child classes."""
        #Make sure temperature doesn't overflow the cloudy table
        if np.max(temp) > 10**8.6:
            temp2 = np.array(temp)
            temp2[np.where(temp2 > 10**8.6)] = 10**8.6
        else:
            temp2 = temp
        return np.float32(self.cloudy_table.ion(elem, ion, den, temp2))

    def _do_interpolation_work(self,pos, vel, elem_den, temp, hh, amumass, line, get_tau):
        """Run the interpolation on some pre-determined arrays, spat out by _read_particle_data"""
        #Factor of 10^-8 converts line width (lambda_X) from Angstrom to cm
        return _Particle_Interpolate(get_tau*1, self.nbins, self.box, self.velfac, self.atime, line.lambda_X*1e-8, line.gamma_X, line.fosc_X, amumass, pos, vel, elem_den, temp, hh, self.axis, self.cofm)

    def particles_near_lines(self, pos, hh,axis=None, cofm=None):
        """Filter a particle list, returning an index list of those near sightlines"""
        if axis == None:
            axis = self.axis
        if cofm == None:
            cofm = self.cofm
        ind = _near_lines(self.box, pos, hh, axis, cofm)
        return ind

    def get_mass_frac(self,elem, data, ind):
        """Get the mass fraction of a given species from a snapshot.
        Arguments:
            elem = name of element
            data = pointer to hdf5 array containing baryons
            ind = index of particles we care about
        Returns mass_frac - mass fraction of this ion
        """
        if elem == "Z":
            mass_frac = np.array(data["GFM_Metallicity"],dtype=np.float32)
        else:
            nelem = self.species.index(elem)
            #Get metallicity of this metal species
            try:
                mass_frac = np.array(data["GFM_Metals"][:,nelem],dtype=np.float32)
            except KeyError:
                #If GFM_Metals is not defined, fall back to primordial abundances
                metal_abund = np.array([0.76, 0.24],dtype=np.float32)
                mass_frac = metal_abund[nelem]
        mass_frac = mass_frac[ind]
        #Deal with floating point roundoff - mass_frac will sometimes be negative
        mass_frac[np.where(mass_frac <= 0)] = 0
        return mass_frac

    def get_particle_number(self, elem, ion, res=8):
        """
        Get the number of particles that contributed significantly
        to the highest column density region (of width proportional
        to the desired spectral resolution).

        This is defined as the number of particles overlapping
        the region with a mass in the ion at least minmass of the
        most massive such particle.

        If the parameter returned is too low, the species is likely unresolved.

        Parameters:
            elem, ion - species to check
            res - spectral resolution, width of region
            minmass - Minimum mass cutoff of particles considered.
        """
        num_important = np.zeros_like(self.axis)
        den = self.get_density(elem, ion)
        for fn in self.files:
            ff = h5py.File(fn, "r")
            data = ff["PartType0"]
            pos = np.array(data["Coordinates"],dtype=np.float32)
            hh = hsml.get_smooth_length(data)
            ff.close()
            #Find particles we care about
            ind = self.particles_near_lines(pos, hh,self.axis,self.cofm)
            pos = pos[ind,:]
            hh = hh[ind]
            del ind
            #For each spectrum find only those particles near the most massive region
            for spec in xrange(self.NumLos):
                #Particles near this spectrum
                ind = self.particles_near_lines(pos, hh, np.array([self.axis[spec],]), np.array([self.cofm[spec,:],]))
                #Largest col. den region
                if np.size(ind) == 0:
                    continue
                maxx = np.where(np.max(den[spec,:])==den[spec,:])[0][0]
                #Region resolution wide around this zone
                region = self.box/self.nbins*np.array(( maxx-res/(2.*self.dvbin), maxx+res/(2.*self.dvbin) ))
                #Need pos. along axis in this region
                axpos = pos[ind,self.axis[spec]-1]
                ind2 = np.where(np.logical_and(axpos - hh[ind] < region[1] , axpos + hh[ind] > region[0]))
                if np.size(ind2) == 0:
                    continue
                num_important[spec]+=np.size(ind2)
        return num_important


    def replace_not_DLA(self, thresh=10**20.3):
        """
        Replace those sightlines which do not contain a DLA with new sightlines, until all sightlines contain a DLA.
        Must implement get_cofm for this to work
        """
        #Declare variables
        found = 0
        wanted = self.NumLos
        cofm_DLA = np.empty_like(self.cofm)
        #Filter
        col_den = self.compute_spectra("H",1,1215,False)
        ind = self.filter_DLA(col_den, thresh)
        H1_DLA = np.empty_like(col_den)
        #Update saves
        top = np.min([wanted, found+np.size(ind)])
        cofm_DLA[found:top] = self.cofm[ind][:top,:]
        H1_DLA[found:top] = col_den[ind][:top,:]
        found += np.size(ind)
        self.discarded = wanted-np.size(ind)
        print "Discarded: ",self.discarded
        while found < wanted:
            #Get a bunch of new spectra
            self.cofm = self.get_cofm()
            col_den = self.compute_spectra("H",1,1215,False)
            ind = self.filter_DLA(col_den, thresh)
            #Update saves
            top = np.min([wanted, found+np.size(ind)])
            cofm_DLA[found:top] = self.cofm[ind][:top-found,:]
            H1_DLA[found:top] = col_den[ind][:top-found,:]
            found += np.size(ind)
            self.discarded += wanted-np.size(ind)
            print "Discarded: ",self.discarded

        #Copy back
        self.cofm=cofm_DLA
        self.colden[("H",1)]=H1_DLA
        #Finalise the cofm array
        self.cofm_final = True

    def get_cofm(self, num = None):
        """Find a bunch more sightlines: should be overridden by child classes"""
        raise NotImplementedError

    def filter_DLA(self, col_den, thresh=10**20.3):
        """Find sightlines with a DLA"""
        #DLAs are huge objects in redshift space (several 10s of A wide), so we want to
        #sum the column densities over the entire spectrum.
        ind = np.where(np.sum(col_den,axis=1) > thresh)
        return ind

    def find_absorber_width(self, elem, ion, chunk = 20, minwidth=None):
        """
           Find the region in velocity space considered to be an absorber for each spectrum.
           This is defined to be the maximum of 1000 km/s and the region over which there is "significant"
           absorption in the strongest line for this ion, where strongest is the line with the largest
           cross-section, ie, greatest lambda * fosc.
           elem, ion - ion to look at

           The threshold above which absorption is considered significant is
           tau = - sigma log(1-F(-3)/SNR). [where F(x) = 1-exp(-x)]
           Probability of a spurious detection from Gaussian noise is N(sigma)**chunk.
           We use N(4) = 1-3e-5, so there is a 0.9994 probability the detection is real for chunk =20.
           If SNR = 0, use 0.2 (for an assumed SNR of 20).

           Returns the low and high indices of absorption, and the offset for the maximal absorption.
        """
        if minwidth == None:
            minwidth = self.minwidth
        try:
            return self.absorber_width[(elem, ion, minwidth)]
        except KeyError:
            pass
        if self.snr > 0:
            thresh = - 4. * np.log(1-(1-np.exp(-3))/self.snr)
        else:
            thresh = 0.2
        lines = self.lines[(elem,ion)]
        strength = [ll.fosc_X*ll.lambda_X for ll in lines.values()]
        ind = np.where(strength == np.max(strength))[0][0]
        #Lines are indexed by wavelength
        strlam = int(lines.values()[ind].lambda_X)
        #Absorption in a strong line: eg, SiII1260.
        strong = self.get_tau(elem, ion, strlam)
        (offset, roll) = self._get_rolled_spectra(strong)
        #Minimum
        if minwidth > 0 and minwidth < self.nbins/2:
            low  = int(self.nbins/2-minwidth/self.dvbin)*np.ones(self.NumLos, dtype=np.int)
            high = int(self.nbins/2+minwidth/self.dvbin)*np.ones(self.NumLos, dtype=np.int)
        else:
            low = np.zeros(self.NumLos, dtype=np.int)
            high = self.nbins*np.ones(self.NumLos, dtype=np.int)
        for ii in xrange(self.NumLos):
            #First expand the search area in case there is absorption at the edges.
            for i in xrange(low[ii],0,-chunk):
                if not np.any(roll[ii,i:(i+chunk)] > thresh):
                    low[ii] = i
                    break
            #Where is there no absorption rightwards of the peak?
            for i in xrange(high[ii],self.nbins,chunk):
                if not np.any(roll[ii,i:(i+chunk)] > thresh):
                    high[ii] = i+chunk
                    break
            #Shrink to width which has some absorption
            ind = np.where(roll[ii][low[ii]:high[ii]] > thresh)[0]
            if np.size(ind) != 0:
                low[ii] = np.max((ind[0]+low[ii]-chunk,0))
                high[ii] = np.min((ind[-1]+low[ii]+chunk,self.nbins))
        self.absorber_width[(elem, ion, minwidth)] = (low, high, offset)
        return (low, high, offset)

    def _eq_width_from_colden(self, col_den, elem = "H", ion = 1, line = 1215):
        """Find the equivalent width of the line from the column density,
           assuming we are in the damping wing regime. Default line is Lyman-alpha.
           Returns width in km/s.
        """
        #line properties
        line = self.lines[(elem,ion)][line]
        #Convert wavelength to cm
        lambdacgs = line.lambda_X*1e-8
        #Tompson cross-section for the electron
        sigma_t = 6.652458558e-25
        #Line cross-section
        sigma_a = np.sqrt(3*math.pi*sigma_t/8.)*lambdacgs*line.fosc_X
        #In the damping wing, W ~ sqrt(N).
        #Coefficients come from setting tau = 1, and expanding the Voigt function used
        #in Tepper-Garcia 2006 where exp(-x^2) ~ 0 (ie, far from the line center)
        #then simplifying a bit
        width = np.sqrt(line.gamma_X*lambdacgs*self.light*col_den*sigma_a)/math.pi
        #Convert from cm/s to km/s
        return width/1e5

    def get_metallicity(self):
        """Return the metallicity, as M/H"""
        MM = self.get_density("Z",-1)
        HH = self.get_density("H",-1)
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
        result =  self._interpolate_single_file(self.files[0], elem, ion, ll, get_tau)
        #Do remaining files
        for fn in self.files[1:]:
            tresult =  self._interpolate_single_file(fn, elem, ion, ll, get_tau)
            #Add new file
            result += tresult
            del tresult
        return result

    def get_observer_tau(self, elem, ion, number=-1, force_recompute=False):
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
            nlines = len(self.lines[(elem,ion)])
            tau = np.empty([nlines, self.NumLos,self.nbins])
            pos = {}
            vel = {}
            elem_den = {}
            temp = {}
            hh = {}
            amumass = {}
            for ff in self.files:
                (pos[ff], vel[ff], elem_den[ff], temp[ff], hh[ff], amumass[ff]) = self._read_particle_data(ff, elem, ion,True)

            for ll in xrange(nlines):
                line = (self.lines[(elem,ion)].values())[ll]
                for ff in self.files:
                    if amumass[ff] != False:
                        tau_loc = self._do_interpolation_work(pos[ff], vel[ff], elem_den[ff], temp[ff], hh[ff], amumass[ff], line, True)
                        tau[ll,:,:] += tau_loc
                        del tau_loc
            #Maximum tau in each spectra with each line,
            #after convolving with a Gaussian for instrumental broadening.
            maxtaus = np.max(self.res_corr(tau,self.spec_res), axis=-1)
            #Array for line indices
            ntau = np.empty([self.NumLos, self.nbins])
            #Use the maximum unsaturated optical depth
            for ii in xrange(self.NumLos):
                # we want unsaturated lines, defined as those with tau < 3
                #which is the maximum tau in the sample of Neeleman 2013
                #Also use lines with some absorption: tau > 0.1, roughly twice noise level.
                ind = np.where(np.logical_and(maxtaus[:,ii] < 3, maxtaus[:,ii] > 0.1))
                if np.size(ind) > 0:
                    line = np.where(maxtaus[:,ii] == np.max(maxtaus[ind,ii]))
                else:
                    #We have no lines in the desired region: here use something slightly saturated.
                    #In reality the observers will use a different ion
                    ind2 = np.where(maxtaus[:,ii] > 0.1)
                    if np.size(ind2) > 0:
                        line = np.where(maxtaus[:,ii] == np.min(maxtaus[ind2,ii]))
                    else:
                        #We have no observable lines: this spectra are metal-poor
                        #and will be filtered anyway.
                        line = np.where(maxtaus[:,ii] == np.max(maxtaus[:,ii]))
                if np.size(line) > 1:
                    line = (line[0][0],)
                ntau[ii,:] = tau[line,ii,:]
            self.tau_obs[(elem, ion)] = ntau
        if number >= 0:
            ntau = ntau[number,:]
        # Convolve lines by a Gaussian filter of the resolution of the spectrograph.
        ntau = self.res_corr(ntau, self.spec_res)
        #Add noise
        if self.snr > 0:
            ntau = self.add_noise(self.snr, ntau, number)
        return ntau

    def res_corr(self, flux, fwhm=8):
        """
           Real spectrographs have finite spectral resolution.
           Correct for this by smoothing the spectrum (the flux) by convolving with a Gaussian.
           The input spectrum is assumed to have infinite resolution, since we have used a spline
           to interpolate it first and/or we are converged.
           Strictly speaking we should rebin the spectrum after to have the same resolution
           as the observed pixel, but as long as the pixels are smaller than the FWHM of the
           spectrograph (which is the case as long as the observer is smart) we will be fine.
           args:
               fwhm - FWHM of the spectrograph in km/s
        """
        # Convert FWHM input to internal units
        res = fwhm/self.dvbin
        #FWHM of a Gaussian is 2 \sqrt(2 ln 2) sigma
        sigma = res/(2*np.sqrt(2*np.log(2)))
        oflux = gaussian_filter1d(flux, sigma, axis=-1)
        return oflux


    def get_filt(self, elem, ion, thresh = 100):
        """
        Get an index list to exclude spectra where the ion is too small, usually the result of
        unresolved star formation.

        thresh - observable density threshold
        """
        #Remember this is not in log.
        met = np.max(self.get_density(elem, ion), axis=1)
        phys = self.dvbin/self.velfac*self.rscale
        ind = np.where(met > thresh/phys)
        return ind

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
                tau_l = np.roll(tau[ll,:],offset[ll])[low[ll]:high[ll]]
                tot_tau = np.sum(tau_l)
                (nnlow, nnhigh) = self._vel_width_bound(tau_l, tot_tau)
                vel_width[ll] = self.dvbin*(nnhigh-nnlow)
            #Return the width
            self.vel_widths[(elem, ion)] = vel_width
            return self.vel_widths[(elem, ion)]

    def _get_rolled_spectra(self,tau):
        """
        Cycle the array tau so that the peak is at the middle.
        """
        tau_out = np.zeros(np.shape(tau))
        roll = np.zeros(np.shape(tau[:,0]), dtype=int)
        for ll in xrange(np.shape(tau)[0]):
            #Deal with periodicity by making sure the deepest point is in the middle
            tau_l = tau[ll,:]
            max_t = np.max(tau_l)
            if max_t == 0:
                continue
            ind_m = np.where(tau_l == max_t)[0][0]
            tau_out[ll] = np.roll(tau_l, np.size(tau_l)/2- ind_m)
            roll[ll] = np.size(tau_l)/2 - ind_m
        return (roll, tau_out)

    def _vel_width_bound(self, tau, tot_tau):
        """Find the 0.05 and 0.95 bounds of the integrated optical depth"""
        cum_tau = np.cumsum(tau)
        #Use spline interpolation to find the edge of the bins.
        tdiff = cum_tau - 0.95*tot_tau
        high = np.where(tdiff >= 0)[0][0]
        tdiff = cum_tau - 0.05*tot_tau
        low = np.where(tdiff >= 0)[0][0]
        return (low, high)

    def _vel_median(self, tau, tot_tau):
        """Find the median point of the integrated optical depth"""
        cum_tau = np.cumsum(tau)
        #Use spline interpolation to find the edge of the bins.
        tdiff = cum_tau - 0.5*tot_tau
        high = np.where(tdiff >= 0)[0][0]
        return high

    def _check_actual_root(self,roots,array):
        """
        Sometimes spl.roots returns multiple roots.
        Ideally we would use a linear spline to avoid this
        but root-finding is not supported for that (?!)
        So check by hand that each of these solutions is an actual root,
        ie that the sign of the array changes there.
        """
        for hh in roots:
            if array[np.floor(hh)]*array[np.ceil(hh)] < 0:
                return hh
        raise ValueError("No root found")

    def vel_mean_median(self, elem, ion):
        """Find the difference between the mean velocity and the median velocity.
           The mean velocity is the point halfway across the extent of the velocity width.
           The median velocity is v(tau = tot_tau /2)
           """
        tau = self.get_observer_tau(elem, ion)
        (low, high, offset) = self.find_absorber_width(elem, ion)
        mean_median = np.zeros(np.shape(tau)[0])
        #Deal with periodicity by making sure the deepest point is in the middle
        for ll in xrange(np.shape(tau)[0]):
            tau_l = np.roll(tau[ll,:],offset[ll])[low[ll]:high[ll]]
            tot_tau = np.sum(tau_l)
            (nnlow, nnhigh) = self._vel_width_bound(tau_l, tot_tau)
            vel_median = self._vel_median(tau_l,tot_tau)
            vmean = (nnlow+nnhigh)/2.
            mean_median[ll] = np.abs(vmean - vel_median)/((nnhigh-nnlow)/2.)
        #Return the width
        return mean_median

    def vel_peak(self, elem, ion):
        """
           Find the f_peak statistic for spectra in an ion.
           f_peak = (vel_peak - vel_mean) / (v_90/2)
        """
        tau = self.get_observer_tau(elem, ion)
        (low, high, offset) = self.find_absorber_width(elem, ion)
        peak = np.zeros(np.shape(tau)[0])
        #Deal with periodicity by making sure the deepest point is in the middle
        for ll in xrange(np.shape(tau)[0]):
            tau_l = np.roll(tau[ll,:],offset[ll])[low[ll]:high[ll]]
            peak[ll] = self._vel_peak_tau(tau_l)
        #Return the width
        return peak

    def _vel_peak_tau(self,tau_l):
        """Helper function for a single spectrum to compute v_peak"""
        tot_tau = np.sum(tau_l)
        (nnlow, nnhigh) = self._vel_width_bound(tau_l, tot_tau)
        vmax = np.where(tau_l == np.max(tau_l))[0][0]
        vmean = (nnlow+nnhigh)/2.
        peak = np.abs(vmax - vmean)/((nnhigh-nnlow)/2.)
        return peak

    def vel_2nd_peak(self, tau):
        """
           Find the difference between the 2nd highest peak optical depth and the mean velocity,
           divided by the velocity width.
           To count as a peak, the 2nd peak must be > 1/3 the peak value,
           and must have a minimum between it and the peak, separated by at least "3-sigma".
           Since these spectra are noiseless, I interpret this as 5%.

           If there is no 2nd peak, return the mean minus the main peak
        """
        #  Size of a single velocity bin
        tot_tau = np.sum(tau,axis = 1)
        mean_median = np.zeros(np.shape(tot_tau))
        tau = self._get_rolled_spectra(tau)
        for ll in np.arange(0, np.shape(tau)[0]):
            (low, high) = self._vel_width_bound(tau[ll,:], tot_tau[ll])
            vmean = low+(high-low)/2.
            #Find second peak
            tt = tau[ll,:][low:high]
            #derivative
            ttd = np.diff(tt)
            x = np.arange(np.size(tt))
            splp = UnivariateSpline(x, ttd, s=0)
            turn = splp.roots()
            spl = UnivariateSpline(x, tt, s=0)
            vals = spl(turn)
            #The peak
            maxpeak = np.where(vals == np.max(vals))
            #The second-highest turning point
            secpeak = np.where(vals == np.max(vals[np.where(vals < vals[maxpeak])]))
            #Is this peak > 1/3 the peak value
            if vals[secpeak] < vals[maxpeak]/3.:
                mean_median[ll] = np.abs(maxpeak+low - vmean)/((high-low)*0.5)
                continue
            #Compute the sign
            sign = -1
            if secpeak < maxpeak and vmean < secpeak:
                sign = 1
            if secpeak > maxpeak and vmean > secpeak:
                sign = 1
            #Find a minimum
            if secpeak < maxpeak:
                minn = np.where((turn < maxpeak)*(turn > secpeak))
            else:
                minn = np.where((turn > maxpeak)*(turn < secpeak))
            #Is the minimum sufficiently deep (and a minimum)
            if np.size(minn == 0) or np.min(vals[minn]) > vals[secpeak]*0.95:
                mean_median[ll] = np.abs(maxpeak+low - vmean)/((high-low)*0.5)
                continue
            mean_median[ll] = sign*np.abs(secpeak+low - vmean)/((high-low)*0.5)
        #Return the width
        return mean_median

    def equivalent_width(self, elem, ion, line):
        """Calculate the equivalent width of a line in Angstroms"""
        tau = self.get_tau(elem, ion, line)
        #1 bin in wavelength: δλ =  λ . v / c
        #λ here is the rest wavelength of the line.
        #speed of light in km /s
        light = self.light / 1e5
        #Line data
        line = self.lines[(elem,ion)][line]
        #lambda in Angstroms, dvbin in km/s,
        #so dl is in Angstrom
        dl = self.dvbin / light * line.lambda_X
        eq_width = cumtrapz(1-np.exp(-tau),dx=dl, axis=1)[:,-1]
        #Don't need to divide by 1+z as lambda_X is already rest wavelength
        return eq_width


    def vel_width_hist(self, elem, ion, dv=0.1):
        """
        Compute a histogram of the velocity widths of our spectra, with the purpose of
        comparing to the data of Prochaska 2008.

        Note this does not match Pontzen 2008, who multiply by the DLA fraction (0.065) obtained from the cddf.

        So we have f(N) = d n/ dv
        and n(N) = number of absorbers per sightline in this velocity bin.
        Note f(N) has dimensions of s/km, because v has units of km/s.

        Parameters:
            elem - element to use
            line - line to use (the components of this line must be pre-computed and stored in self.metals)
            dv - bin spacing
            met_cut - Discard spectra whose maximal metal column density is below this level.
                      Removes unobservable systems.

        Returns:
            (v, f_table) - v (binned in log) and corresponding f(N)
        """
        return self._vel_stat_hist(elem, ion, dv, self.vel_width, log=True)

    def f_meanmedian_hist(self, elem, ion, dv=0.1):
        """
        Compute a histogram of the mean median statistic of our spectra, the difference in
        units of the velocity width between the mean velocity and median velocity of
        the absorber.

        For arguments see vel_width_hist.
        """
        return self._vel_stat_hist(elem, ion, dv, self.vel_mean_median, log=False)

    def f_peak_hist(self, elem, ion, dv=0.1):
        """
        Compute a histogram of the peak statistic of our spectra, the difference in
        units of the velocity width between the largest peak velocity and the mean velocity of
        the absorber.

        For arguments see vel_width_hist.
        """
        return self._vel_stat_hist(elem, ion, dv, self.vel_peak, log=False)

    def eq_width_hist(self, elem, ion, line, dv=0.05, eq_cut = 0.02):
        """
        Compute a histogram of the equivalent width distribution of our spectra, with the purpose of
        comparing to the data of Neeleman 2013.

        Returns:
            (v, f_table) - v (binned in log) and corresponding f(N)
        """
        print "For ",self.lines[(elem,ion)][line].lambda_X," Angstrom"
        eq_width = self.equivalent_width(elem, ion, line)
        #Filter small eq. widths as they would not be detected
        ind = np.where(eq_width > eq_cut)
        eq_width = eq_width[ind]
        v_table = np.arange(np.log10(np.min(eq_width)), np.log10(np.max(eq_width)), dv)
        vbin = np.array([(v_table[i]+v_table[i+1])/2. for i in range(0,np.size(v_table)-1)])
        eqws = np.histogram(np.log10(eq_width),v_table, density=True)[0]
        return (vbin, eqws)

    def _vel_stat_hist(self, elem, ion, dv, func, log=True, filt=True):
        """
           Internal function that finds the histogram in velocity space of
           the values of a statistic for a particular ion.
        """
        #Filter small number of spectra without metals
        vel_width = func(elem, ion)
        if filt:
            ind = self.get_filt(elem, ion)
            vel_width = vel_width[ind]
        if np.size(dv) > 1:
            v_table = dv
        elif log:
            v_table = 10**np.arange(0, np.min((50,np.log10(np.max(vel_width)))), dv)
        else:
            v_table = np.arange(0, 1, dv)
        vbin = np.array([(v_table[i]+v_table[i+1])/2. for i in range(0,np.size(v_table)-1)])
        if log:
            vels = np.histogram(np.log10(vel_width),np.log10(v_table), density=True)[0]
        else:
            vels = np.histogram(vel_width,v_table, density=True)[0]
        return (vbin, vels)

    def mass_hist(self, dm=0.1):
        """
        Compute a histogram of the host halo mass of each DLA spectrum.

        Parameters:
            dm - bin spacing

        Returns:
            (mbins, pdf) - Mass (binned in log) and corresponding PDF.
        """
        (halos, _) = self.find_nearest_halo()
        f_ind = np.where(halos != -1)
        #nlos = np.shape(vel_width)[0]
        #print 'nlos = ',nlos
        virial = self.virial_vel(halos[f_ind])
        m_table = 10**np.arange(np.log10(np.min(virial)+0.1), np.log10(np.max(virial)), dm)
        mbin = np.array([(m_table[i]+m_table[i+1])/2. for i in range(0,np.size(m_table)-1)])
        pdf = np.histogram(np.log10(virial),np.log10(m_table), density=True)[0]
        print "Field DLAs: ",np.size(halos)-np.size(f_ind)
        return (mbin, pdf)

    def virial_vel(self, halos=None, subhalo=False):
        """Get the virial velocities of the selected halos in km/s"""
        if subhalo:
            if halos != None:
                mm = self.sub_sub_mass[halos]
                rr = np.array(self.sub_sub_radii[halos])
            else:
                mm = self.sub_sub_mass
                rr = np.array(self.sub_sub_radii)
        else:
            if halos != None:
                mm = self.sub_mass[halos]
                rr = np.array(self.sub_radii[halos])
            else:
                mm = self.sub_mass
                rr = np.array(self.sub_radii)
        #physical cm from comoving kpc/h
        cminkpch = self.UnitLength_in_cm/self.hubble/(1+self.red)
        # Conversion factor from M_sun/kpc/h to g/cm
        conv = self.UnitMass_in_g/1e10/cminkpch
        #Units: grav is in cm^3 /g/s^-2
        #Define zero radius, zero mass halos as having zero virial velocity.
        rr[np.where(rr == 0)] = 1
        virial = np.sqrt(self.gravcgs*conv*mm/rr)/1e5
        return virial

    def get_col_density(self, elem, ion, force_recompute=False):
        """get the column density in each pixel for a given species"""
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
        """Get the density in each pixel for a given species"""
        colden = self.get_col_density(elem, ion, force_recompute)
        phys = self.dvbin/self.velfac*self.rscale
        return colden/phys

    def get_tau(self, elem, ion,line, number = -1, force_recompute=False):
        """Get the column density in each pixel for a given species"""
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
        tau = self.res_corr(tau, self.spec_res)
        if self.snr > 0:
            tau = self.add_noise(self.snr, tau, number)
        return tau

    def _vel_single_file(self,fn, elem, ion):
        """Get the column density weighted interpolated velocity field for a single file"""
        (pos, vel, elem_den, temp, hh, amumass) = self._read_particle_data(fn, elem, ion,True)
        if amumass == False:
            return np.zeros([np.shape(self.cofm)[0],self.nbins,3],dtype=np.float32)
        else:
            line = self.lines[("H",1)][1215]
            vv =  np.empty([np.shape(self.cofm)[0],self.nbins,3],dtype=np.float32)
            phys = self.dvbin/self.velfac*self.rscale
            for ax in (0,1,2):
                weight = vel[:,ax]*np.sqrt(self.atime)
                vv[:,:,ax] = self._do_interpolation_work(pos, vel, elem_den*weight/phys, temp, hh, amumass, line, False)
            return vv

    def get_velocity(self, elem, ion):
        """Get the column density weighted velocity in each pixel for a given species.
        """
        try:
            self._really_load_array((elem, ion), self.velocity, "velocity")
            return self.velocity[(elem, ion)]
        except KeyError:
            velocity =  self._vel_single_file(self.files[0], elem, ion)
            #Do remaining files
            for fn in self.files[1:]:
                tresult =  self._vel_single_file(fn, elem, ion)
                #Add new file
                velocity += tresult
                del tresult
            den = self.get_density(elem, ion)
            den[np.where(den == 0.)] = 1
            #Broadcasting can't handle this
            for ax in (0,1,2):
                velocity[:,:,ax] /= den
            self.velocity[(elem, ion)] = velocity
            return velocity

    def _temp_single_file(self,fn, elem, ion):
        """Get the density weighted interpolated temperature field for a single file"""
        (pos, vel, elem_den, temp, hh, amumass) = self._read_particle_data(fn, elem, ion,True)
        if amumass == False:
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
            temp =  self._temp_single_file(self.files[0], elem, ion)
            #Do remaining files
            for fn in self.files[1:]:
                tresult =  self._temp_single_file(fn, elem, ion)
                #Add new file
                temp += tresult
                del tresult
            den = self.get_density(elem, ion)
            den[np.where(den == 0.)] = 1
            temp /= den
            self.temp[(elem, ion)] = temp
            return temp

    def column_density_function(self,elem = "H", ion = 1, dlogN=0.2, minN=13, maxN=23.):
        """
        This computes the DLA column density function, which is the number
        of absorbers per sight line with HI column densities in the interval
        [NHI, NHI+dNHI] at the absorption distance X.
        Absorption distance is simply a single simulation box.
        A sightline is assumed to be equivalent to one grid cell.
        That is, there is presumed to be only one halo in along the sightline
        encountering a given halo.

        So we have f(N) = d n_DLA/ dN dX
        and n_DLA(N) = number of absorbers per sightline in this column density bin.
                     1 sightline is defined to be one grid cell.
                     So this is (cells in this bins) / (no. of cells)
        ie, f(N) = n_DLA / ΔN / ΔX
        Note f(N) has dimensions of cm^2, because N has units of cm^-2 and X is dimensionless.

        Parameters:
            dlogN - bin spacing
            minN - minimum log N
            maxN - maximum log N

        Returns:
            (NHI, f_N_table) - N_HI (binned in log) and corresponding f(N)
        """
        NHI_table = 10**np.arange(minN, maxN, dlogN)
        center = np.array([(NHI_table[i]+NHI_table[i+1])/2. for i in range(0,np.size(NHI_table)-1)])
        width =  np.array([NHI_table[i+1]-NHI_table[i] for i in range(0,np.size(NHI_table)-1)])
        dX=self.absorption_distance()/self.nbins
        #Col density of each pixel: makes more sense for lower column densities
        rho = np.ravel(self.get_col_density(elem, ion))
        tot_cells = np.size(rho)
        (tot_f_N, NHI_table) = np.histogram(rho,NHI_table)
        tot_f_N=tot_f_N/(width*dX*tot_cells)
        return (center, tot_f_N)

    def absorption_distance(self):
        """
        Compute X(z), the absorption distance per sightline (dimensionless)
        X(z) = int (1+z)^2 H_0 / H(z) dz
        When dz is small, dz ~ H(z)/c dL, so
        X(z) ~ (1+z)^2 H_0/c dL
        """
        #h * 100 km/s/Mpc in h/s
        h100=3.2407789e-18
        #Units: h/s   s/cm                 kpc/h      cm/kpc
        return h100/self.light*self.box*self.UnitLength_in_cm*(1+self.red)**2

    def rho_crit(self):
        """Get the critical density at z=0 in units of g cm^-3"""
        #H in units of 1/s
        h100=3.2407789e-18*self.hubble
        #G in cm^3 g^-1 s^-2
        grav=6.672e-8
        rho_crit=3*h100**2/(8*math.pi*grav)
        return rho_crit

    def _rho_DLA(self, thresh=10**20.3, elem = "H", ion = 1):
        """Compute rho_DLA, the sum of the mass in DLAs,
           divided by the volume of the spectra in g/cm^3 (comoving).
            ρ_DLA = m_p * avg. column density / (1+z)^2 / length of column
            Note: If we want the neutral gas density rather than the
            neutral hydrogen density, divide by 0.76, the hydrogen mass fraction.
        """
        #Column density of HI in atoms cm^-2 (physical)
        col_den = self.get_col_density(elem, ion)
        if thresh > 0:
            HIden = np.sum(col_den[np.where(col_den > thresh)])/np.size(col_den)
        else:
            HIden = np.mean(col_den)
        HIden *= np.size(col_den)/(np.size(col_den)+1.*self.discarded*self.nbins)
        #Avg. Column density of HI in g cm^-2 (comoving)
        HIden = self.protonmass * HIden/(1+self.red)**2
        #Length of column (each cell) in comoving cm
        length = (self.box*self.UnitLength_in_cm/self.hubble)/self.nbins/(1+self.red)
        #Avg density in g/cm^3 (comoving)
        return HIden/length

    def rho_DLA(self, thresh=10**20.3):
        """Compute rho_DLA, the average density in DLAs. This is almost the same as the average density in HI.
           Units are 10^8 M_sun / Mpc^3 (comoving), like 0811.2003
        """
        #Avg density in g/cm^3 (comoving)
        rho_DLA = self._rho_DLA(thresh)
        # 1e8 M_sun/Mpc^3 in g/cm^3
        conv = 0.01 * self.UnitMass_in_g / self.UnitLength_in_cm**3
        return rho_DLA / conv

    def omega_DLA(self, thresh=10**20.3, elem = "H", ion = 1):
        """Compute Omega_DLA, the sum of the mass in DLAs, divided by the volume of the spectra, divided by the critical density.
            Ω_DLA = m_p * avg. column density / (1+z)^2 / length of column / rho_c
            Note: If we want the neutral gas density rather than the neutral hydrogen density, divide by 0.76,
            the hydrogen mass fraction.
        """
        #Avg density in g/cm^3 (comoving) divided by critical density in g/cm^3
        omega_DLA=self._rho_DLA(thresh, elem, ion)/self.rho_crit()
        return omega_DLA

    def line_density(self, thresh=10**20.3, elem = "H", ion = 1):
        """Compute the line density, the total no. of DLA sightlines divided by the total number of sightlines, multiplied by d L / dX. This is dN/dX = l_DLA(z)
        """
        col_den = self.get_col_density(elem, ion)
        #Average fraction of pixels containing a DLA
        frac = 1.*np.size(col_den[np.where(col_den > thresh)])/np.size(col_den)
        #Divide by abs. distance per sightline
        frac *= np.size(col_den)/(np.size(col_den)+1.*self.discarded*self.nbins)
        return frac/(self.absorption_distance()/self.nbins)

    def get_separated(self, elem="Si", ion = 2, thresh = 1e-1, mindist=0):
        """
        Find spectra with more than a single density peak.
        Threshold is as a percentage of the maximum value.
        mindist is in km/s
        """
        dist = int(mindist/self.dvbin)
        ind = self.get_filt(elem, ion)
        rho = self.get_col_density(elem, ion)[ind]
        H1_den = self.get_col_density("H", 1)[ind]
        seps = np.zeros(np.size(ind[0]), dtype=np.bool)
        lls = 0
        dla = 0
        none = 0
        #deal with periodicity by making sure the deepest point is in the middle
        for ll in xrange(np.size(ind[0])):
            rho_l = rho[ll,:]
            H1_l = H1_den[ll,:]
            lsep = combine_regions(rho_l > thresh*np.max(rho_l), dist)
            seps[ll] = (np.shape(lsep)[0] > 1)
            if seps[ll] == False:
                continue
            m_H1_colden = np.array([np.sum(H1_l[lsep[jj,0]:lsep[jj,1]]) for jj in xrange(np.shape(lsep)[0])])
            #All DLAs
            if np.all(m_H1_colden > 10**(20.3)):
                dla += 1
            #Some LLS
            elif np.all(m_H1_colden > 10**(17.)):
                lls += 1
            else:
                none += 1
        tot = dla + lls + none
        print "Fraction DLA: ",1.*dla/tot," Fraction LLS: ",1.*lls/tot," fraction less: ",1.*none/tot
        return seps

    def get_spectra_proj_pos(self, cofm=None):
        """Get the position of the spectra projected to their origin"""
        if np.mean(self.axis) != self.axis[0] or  self.axis[0] != self.axis[-1]:
            raise ValueError("Not all spectra are along the same axis")
        if cofm == None:
            cofm = self.cofm
        axis = self.axis[0]
        if axis == 1:
            spos = cofm[:,1:]
        if axis == 2:
            spos = np.vstack([cofm[:,0],cofm[:,2]]).T
        if axis == 3:
            spos = cofm[:,:2]
        return spos

    def min_halo_mass(self, minpart = 400):
        """Min resolved halo mass in internal Gadget units (1e10 M_sun)"""
        #This is rho_c in units of h^-1 1e10 M_sun (kpc/h)^-3
        rhom = 2.78e+11* self.OmegaM / 1e10 / (1e3**3)
        #Mass of an SPH particle, in units of 1e10 M_sun, x omega_m/ omega_b.
        try:
            target_mass = self.box**3 * rhom / self.npart[0]
        except AttributeError:
            #Back-compat hack
            target_mass = self.box**3 * rhom / 512.**3
        min_mass = target_mass * minpart
        return min_mass

    def load_halo(self):
        """Load a halo catalogue: note this will return some halos with zero radius.
           These are empty FoF groups and should not be a problem."""
        SolarMass_in_g=1.989e33
        try:
            subs=subfindhdf.SubFindHDF5(self.base, self.num)
            #Get particle center of mass, use group catalogue.
            self.sub_cofm=subs.get_grp("GroupPos")
            #halo masses in M_sun/h: use M_200
            self.sub_mass=subs.get_grp("Group_M_Crit200")*self.UnitMass_in_g/SolarMass_in_g
            #r200 in kpc/h (comoving).
            self.sub_radii = subs.get_grp("Group_R_Crit200")
            self.sub_vel = subs.get_grp("GroupVel")
            self.sub_sub_radii =  subs.get_sub("SubhaloHalfmassRad")
            self.sub_sub_cofm =  subs.get_sub("SubhaloPos")
            self.sub_sub_mass =  subs.get_sub("SubhaloMass")
            self.sub_sub_index = subs.get_sub("SubhaloGrNr")
            self.sub_sub_vel = subs.get_sub("SubhaloVel")
        except IOError:
            pass

    def assign_to_halo(self, zpos, halo_radii, halo_cofm):
        """
        Assign a list of lists of positions to halos, by finding the unique halo
        within whose virial radius each position is.
        """
        dists = []
        halos = []
        #X axis first
        for ii in xrange(len(zpos)):
            proj_pos = np.array(self.cofm[ii,:])
            ax = self.axis[ii]-1
            dists.append([])
            halos.append([])
            for zzp in zpos[ii]:
                proj_pos[ax] = zzp
                #Is this within the virial radius of any halo?
                dd = ne.evaluate("sum((halo_cofm - proj_pos)**2,axis=1)")
                ind = np.where(dd < halo_radii**2)
                #Should not be multiple close halos
                # assert(np.size(ind) < 2)
                #Very rarely, in 2/5000 cases,
                #something hits the edge of more than one halo
                #This is so rare we don't worry about it.
                if np.size(ind) >= 1:
                    halos[ii].append(ind[0][0])
                    dists[ii].append(np.sqrt(dd[ind][0]))
        return (halos, dists)

    def get_contiguous_regions(self, elem="H", ion = 1, thresh = 2e20, relthresh = 1e-3):
        """
        Find the weighted z position of all contiguous DLA-hosting regions in each spectrum.
        Returns a list of lists. Each element in the outer list corresponds to a spectrum.
        Each inner list is the list of weighted z positions of DLA-hosting regions.
        """
        den = self.get_density(elem, ion)
        contig = []
        seps = np.zeros(self.NumLos, dtype=np.bool)
        (roll, colden) = self._get_rolled_spectra(den)
        #deal with periodicity by making sure the deepest point is in the middle
        for ii in xrange(self.NumLos):
            # This is column density, not absorption, so we cannot
            # use the line width to find the peak region.
            # Instead, just search a constant velocity offset from the peak,
            #to roughly match the velocity width calculation.
            lcolden = colden[ii,self.nbins/2-1000*self.dvbin:self.nbins/2+1000*self.dvbin]
            offset = self.nbins/2-1000*self.dvbin - roll[ii]
            # Get first and last indices of separate regions in list
            if np.max(lcolden) > thresh:
                seps = combine_regions(lcolden > thresh)
            else:
                seps = combine_regions(lcolden > relthresh*np.max(lcolden))
            # Find weighted z position for each one
            zposes = []
            for jj in xrange(np.shape(seps)[0]):
                nn = offset+np.arange(self.nbins)[seps[jj,0]:seps[jj,1]]
                llcolden = lcolden[seps[jj,0]:seps[jj,1]]
                zpos = ne.evaluate("sum(llcolden*nn)")
                summ = ne.evaluate("sum(llcolden)")
                #Make sure it refers to a valid position
                zpos = (zpos / summ) % self.nbins
                zpos *= 1.*self.box/self.nbins
                zposes.append(zpos)
            contig.append(zposes)
        return contig

    def find_nearest_halo(self):
        """Find the single most massive halos associated with absorption near a sightline, possibly via a subhalo."""
        (halos, subhalos) = self.find_nearby_halos()
        outhalos = np.zeros(self.NumLos,dtype=int)-1
        for ii in xrange(self.NumLos):
            subhalo_parent = list(self.sub_sub_index[subhalos[ii]])
            both = list(set(subhalo_parent+halos[ii]))
            if len(both) > 0:
                vir_vel = self.virial_vel(both)
                ind = np.where(vir_vel == np.max(vir_vel))
                outhalos[ii] = both[ind[0][0]]
            else:
                outhalos[ii] = -1
        return (outhalos, 0)

    def find_nearby_halos(self):
        """Find halos and subhalos associated with absorption near a sightline"""
        try:
            return (self.spectra_halos, self.spectra_subhalos)
        except AttributeError:
            pass
        zpos = self.get_contiguous_regions(thresh = 1e19, relthresh = 1e-2)
        (halos, _) = self.assign_to_halo(zpos, self.sub_radii, self.sub_cofm)
        (subhalos, _) = self.assign_to_halo(zpos, self.sub_sub_radii, self.sub_sub_cofm)
        #Merge absorption features inside the same halo
        for ii in xrange(self.NumLos):
            halos[ii] = list(set(halos[ii]))
            subhalos[ii] = list(set(subhalos[ii]))
        print "no. halos: ",sum([len(hh) for hh in halos])," mult halos: ",sum([len(hh) > 1 for hh in halos])
        print "no. subhalos: ",sum([len(hh) for hh in subhalos])," mult subhalos: ",sum([len(hh) > 1 for hh in subhalos])
        self.spectra_halos = halos
        self.spectra_subhalos = subhalos
        return (halos, subhalos)

def combine_regions(condition, mindist=0):
    """Combine contiguous regions that are shorter than mindist"""
    reg = contiguous_regions(condition)
    #Find lengths to ignore
    if mindist > 0 and np.shape(reg)[0] > 1:
        newreg = np.array(reg[0,:])
        newreg.shape = (1,2)
        for ii in xrange(1,np.shape(reg)[0]):
            if reg[ii,0] - newreg[-1,1] < mindist:
                #Move the end point of the last segment to that of this one
                newreg[-1,1] = reg[ii,1]
            else:
                #This segment is far from the last one.
                #Add the new segment to the list
                newreg = np.vstack([newreg, reg[ii,:]])
        reg = newreg
    return reg

def contiguous_regions(condition):
    """Finds contiguous True regions of the boolean array "condition". Returns
    a 2D array where the first column is the start index of the region and the
    second column is the end index.
    If mindist != 0, ignores changes shorter than mindist
    """
    # Find the indicies of changes in "condition"
    d = np.diff(condition)
    idx, = d.nonzero()
    # We need to start things after the change in "condition". Therefore,
    # we'll shift the index by 1 to the right.
    idx += 1

    if condition[0]:
        # If the start of condition is True prepend a 0
        idx = np.r_[0, idx]

    if condition[-1]:
        # If the end of condition is True, append the length of the array
        idx = np.r_[idx, condition.size]

    # Reshape the result into two columns
    idx.shape = (-1,2)
    return idx
