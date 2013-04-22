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
import line_data
import h5py
import hdfsim
from scipy.interpolate import UnivariateSpline
import os.path as path
from _spectra_priv import _SPH_Interpolate, _near_lines,_Compute_Absorption

class Spectra:
    """Class to interpolate particle densities along a line of sight and calculate their absorption
        Arguments:
            num - Snapshot number
            base - Name of base directory for snapshot
            cofm - table of los positions, as [n, 3] shape array.
            axis - axis along which to put the sightline
            res (optional) - Spectra pixel resolution in km/s
    """
    def __init__(self,num, base,cofm, axis, res=1., savefile="spectra.hdf5", savedir=None):
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
        if savedir == None:
            savedir=base
        self.savefile = path.join(savedir,"snapdir_"+str(self.num).rjust(3,'0'),savefile)
        #Snapshot data
        try:
            try:
                self.cofm
            except AttributeError:
                self.load_savefile(self.savefile)
        except (IOError, KeyError):
            print "Reloading from snapshot"
            self.cofm = cofm
            self.axis = np.array(axis, dtype = np.int32)
            self.files = hdfsim.get_all_files(num, base)
            ff = h5py.File(self.files[0], "r")
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
            #Empty dictionary to add results to
            self.metals = {}

        # Conversion factors from internal units
        rscale = (self.KPC*self.atime)/self.hubble    # convert length to m
        mscale = (1.0e10*self.SOLAR_MASS)/self.hubble   # convert mass to kg
        self.dscale = mscale / rscale **3 # Convert density to kg / m^3
        #  Calculate the length scales to be used in the box
        self.Hz = 100.0*self.hubble * np.sqrt(self.OmegaM/self.atime**3 + self.OmegaLambda)
        self.vmax = self.box * self.Hz * rscale/ MPC # box size (kms^-1)
        try:
            self.dzgrid   = self.box * rscale / (1.*self.nbins) # bin size m
            self.dvbin = self.dzgrid * self.Hz / MPC # velocity bin size (kms^-1)
        except AttributeError:
            #This will occur if we are not reloading from a snapshot
            self.dvbin = res # velocity bin size (kms^-1)
            self.dzgrid = self.dvbin * MPC / self.Hz #bin size m
            self.nbins = int(self.box * rscale / self.dzgrid) #Number of bins to achieve the required resolution
        #Species we can use
        self.species = ['H', 'He', 'C', 'N', 'O', 'Ne', 'Mg', 'Si', 'Fe']
        #Generate cloudy tables
        self.cloudy_table = convert_cloudy.CloudyTable(self.red)
        #Line data
        self.lines = line_data.LineData()
        print np.size(self.axis), " sightlines. resolution: ", self.dvbin, " z=", self.red

    def save_file(self):
        """
        Saves spectra to a file, because they are slow to generate.
        File is by default to be $snap_dir/snapdir_$snapnum/spectra.hdf5.
        """
        try:
            f=h5py.File(self.savefile,'a')
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
        grp = f.create_group("spectra")
        grp["cofm"]=self.cofm
        grp["axis"]=self.axis
        grp_grid = f.create_group("metals")
        for (key, value) in self.metals.iteritems():
            try:
                #Select by the element
                gg = grp_grid[key[0]]
            except KeyError:
                grp_grid.create_group(key[0])
                gg = grp_grid[key[0]]
            gg.create_dataset(str(key[1]),data=value)
        f.close()

    def load_savefile(self,savefile=None):
        """Load data from a file"""
        #Name of savefile
        f=h5py.File(savefile,'r')
        grid_file=f["Header"]
        self.red=grid_file.attrs["redshift"]
        self.atime = 1./(1+self.red)
        self.OmegaM=grid_file.attrs["omegam"]
        self.nbins=grid_file.attrs["nbins"]
        self.omegab=grid_file.attrs["omegab"]
        self.OmegaLambda=grid_file.attrs["omegal"]
        self.hubble=grid_file.attrs["hubble"]
        self.box=grid_file.attrs["box"]
        self.metals = {}
        grp = f["metals"]
        for elem in grp.keys():
            for ion in grp[elem].keys():
                self.metals[(elem, int(ion))] = np.array(grp[elem][ion])
        grp=f["spectra"]
        self.cofm = np.array(grp["cofm"])
        self.axis = np.array(grp["axis"])
        f.close()

    def SPH_Interpolate_metals(self, elem, ion, get_rho_H=False, ind=None):
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
        (rho_H, rho_metal, vel_metal, temp_metal) =  self._interpolate_single_file(self.files[0], elem, ion, get_rho_H, ind)
        #Do remaining files
        for fn in self.files[1:]:
            (trho_H, trho_metal, tvel_metal, ttemp_metal) =  self._interpolate_single_file(fn, elem, ion, get_rho_H, ind)
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


    def _interpolate_single_file(self,fn, elem, ion, rho_H, h_ind=None):
        """Read arrays and perform interpolation for a single file"""
        ff = h5py.File(fn, "r")
        data = ff["PartType0"]
        pos = np.array(data["Coordinates"],dtype=np.float32)
        hh = hsml.get_smooth_length(data)
        #Find particles we care about
        ind = self.particles_near_lines(pos, hh)
        pos = pos[ind,:]
        hh = hh[ind]
        #Get the rest of the arrays: reducing them each time to have a smaller memory footprint
        vel = np.array(data["Velocities"],dtype=np.float32)
        vel = vel[ind,:]
        mass = np.array(data["Masses"],dtype=np.float32)
        mass = mass[ind]
        u = np.array(data["InternalEnergy"],dtype=np.float32)
        u = u[ind]
        ne = np.array(data["ElectronAbundance"],dtype=np.float32)
        ne = ne[ind]
        ff.close()
        metal_in = self.get_mass_frac(fn, elem, ion, ind)
        #for xx in [pos, vel, mass, u, ne, hh]:
        #    if np.size(np.where(np.isnan(xx))[0]) > 0:
        #        raise ValueError
        #Get rid of ind so we have some memory for the interpolator
        del ind
        if h_ind != None:
            cofm = self.cofm[h_ind]
            axis = self.axis[h_ind]
        else:
            cofm = self.cofm
            axis = self.axis
        out =  _SPH_Interpolate(rho_H*1,self.nbins, self.box, pos, vel, mass, u, ne, metal_in, hh, axis, cofm)
        if not rho_H:
            out = (None,)+out
        return out

    def particles_near_lines(self, pos, hh):
        """Filter a particle list, returning an index list of those near sightlines"""
        ind = _near_lines(self.box, pos, hh, self.axis, self.cofm)
        return ind

    def get_mass_frac(self,fn,elem, ion,ind=None):
        """Get the mass fraction of a given species from a snapshot.
        Arguments:
            fn = file to read
            elem = name of element
            ion = ion species
            ind = optional pre-computed index list of particles we care about
        Returns mass_frac - mass fraction of this ion
        """
        nelem = self.species.index(elem)
        ff = h5py.File(fn,"r")
        data = ff["PartType0"]
        if ind==None:
            pos = np.array(data["Coordinates"],dtype=np.float32)
            hh = np.array(hsml.get_smooth_length(data),dtype=np.float32)
            #Find particles we care about
            ind = self.particles_near_lines(pos, hh)

        #Get metallicity of this metal species
        try:
            mass_frac = np.array(data["GFM_Metals"][:,nelem],dtype=np.float32)
            mass_frac = mass_frac[ind]
            #Deal with floating point roundoff - mass_frac will sometimes be negative
            #10^-30 is Cloudy's definition of zero.
            mass_frac[np.where(mass_frac < 1e-30)] = 1e-30
        except KeyError:
            #If GFM_Metals is not defined, fall back to a guess.
            #Some default abundances. H and He are primordial, the rest are Milky Way as given by wikipedia
            metal_abund = np.array([0.76, 0.24, 4.6e-3, 9.6e-4, 1.04e-2, 1.34e-3, 5.8e-4, 6.5e-4, 1.09e-3],dtype=np.float32)
            mass_frac = metal_abund[nelem]*np.ones(np.shape(data["Density"]),dtype=np.float32)

        #In kg/m^3
        den = np.array(data["Density"], dtype = np.float32)*self.dscale
        #In (hydrogen) atoms / cm^3
        den /= (self.PROTONMASS*100**3)
        den = den[ind]
        #Get density of this ion - we need to weight T and v by ion abundance
        #Cloudy density in physical H atoms / cm^3
        #Special case H1:
        if elem == 'H':
            if ion != 1:
                raise ValueError
            # Hydrogen mass frac in the data array
            mass_frac *= np.array(data["NeutralHydrogenAbundance"],dtype=np.float32)[ind]
        elif ion != -1:
            mass_frac *= self.cloudy_table.ion(elem, ion, mass_frac, den)
        ff.close()
        return mass_frac


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
            profile[ind] = T1[ind] - aa_H1[ind]/np.sqrt(math.pi)/T0[ind]*(T1[ind]**2*(4.0*T0[ind]**2 + 7.0*T0[ind] + 4.0 + T2[ind]) - T2[ind] -1.0)
            tau[i] = np.sum(A_H1  * rho  * profile /(mass*self.PROTONMASS*bb))

        return tau

    def get_tau(self, elem, ion, ll=0):
        """Get the optical depth for a particular element out of:
           (He, C, N, O, Ne, Mg, Si, Fe)
           and some ion number
        """
        try:
            return self.metals[(elem, ion)][3]
        except KeyError:
            #generate metal and hydrogen spectral densities
            #Indexing is: rho_metals [ NSPECTRA, NBIN ]
            [rho, vel, temp] = self.SPH_Interpolate_metals(elem, ion)

            #Compute tau for this metal ion
            (nlos, nbins) = np.shape(rho)
            tau = np.array([self.compute_absorption(elem, ion, ll, rho[n,:], vel[n,:], temp[n,:]) for n in xrange(0, nlos)])
            self.metals[(elem, ion)] = [rho, vel, temp, tau]
            return tau
        except IndexError:
            #This occurs when we have calculated rho, vel and T, but not tau
            [rho, vel, temp] = self.metals[(elem, ion)]
            #Compute tau for this metal ion
            (nlos, nbins) = np.shape(rho)
            tau = np.array([self.compute_absorption(elem, ion, ll, rho[n,:], vel[n,:], temp[n,:]) for n in xrange(0, nlos)])
            self.metals[(elem, ion)] = [rho, vel, temp, tau]
            return tau

    def get_filt(self, elem, line, HI_cut = 10**17, met_cut = 1e13):
        """
        Get an index list of spectra with a DLA in them, and metal column density above a certain value
        """
        #Remember this is not in log...
        if HI_cut == None:
            HI_cut = 0
        if met_cut == None:
            met_cut = 0
        rho = self.get_col_density(elem,line)
        rho_H = self.get_col_density("H",1)
        ind = np.where(np.logical_and(np.max(rho,axis=1) > met_cut, np.max(rho_H,axis=1) > HI_cut))
        return ind

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
            max_t = np.max(tau_l)
            if max_t == 0:
                vel_width[ll] = 0
                continue
            ind_m = np.where(tau_l == max_t)[0][0]
            tau_l = np.roll(tau_l, np.size(tau_l)/2- ind_m)
            (low, high) = self._vel_width_bound(tau_l, tot_tau[ll])
            vel_width[ll] = self.dvbin*(high-low)
        #Return the width
        return vel_width

    def _vel_width_bound(self, tau, tot_tau):
        """Find the 0.05 and 0.95 bounds of the integrated optical depth"""
        cum_tau = np.cumsum(tau)
        #Use spline interpolation to find the edge of the bins.
        tdiff = cum_tau - 0.95*tot_tau
        x = np.arange(0,np.size(cum_tau))
        spl = UnivariateSpline(x, tdiff, s=0)
        high = spl.roots()
        tdiff = cum_tau - 0.05*tot_tau
        spl = UnivariateSpline(np.arange(0,np.size(cum_tau)), tdiff, s=0)
        low = spl.roots()
        if np.size(low) == 0:
            low = 0
        if np.size(high) == 0:
            high = np.size(cum_tau)-1
        if np.size(low) > 1:
            low = low[0]
        if np.size(high) > 1:
            high = high[0]

        return (low, high)

    def delta(self, rho):
        """Get a density in units of delta = ρ/bar{ρ} -1.
        Supplied density should be in physical kg/m^3."""
        #Gravitational constant in SI
        GRAVITY = 6.67428e-11
        # 100kms^-1Mpc^-1 in SI
        H0 = 1.0e5/(self.KPC*1e3)
        # Critical matter/energy density at z = 0
        rhoc = 3.0 * (H0*self.hubble)*(H0*self.hubble) / (8.0 * math.pi * GRAVITY)
        #Primordial hydrogen mass-fraction
        XH = 0.76
        #Mean hydrogen mass density of the Universe in kg/m^3
        critH = (rhoc * self.omegab * XH) / self.atime**3

        # H density normalised by mean
        return rho/critH

    def vel_width_hist(self, elem, line, dv=0.1, HI_cut = None, met_cut = 1e13, unres = 5):
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
            HI_cut - Prochaska used a subsample of spectra containing a DLA.
                     If this value is not None, consider only HI column densities above this threshold.
                     If the spectra are taken within the halo virial radius, this does not make much of a difference.
            met_cut - Discard spectra whose maximal metal column density is below this level.
                      Removes unobservable systems.
            unres - Remove systems with velocity widths below this value, where they are affected
                    by the pixel size of the spectra.

        Returns:
            (v, f_table) - v (binned in log) and corresponding f(N)
        """
        tau = self.metals[(elem, line)][3]

        vel_width = self.vel_width(tau[self.get_filt(elem, line, HI_cut, met_cut)])
        if unres != None:
            ind = np.where(vel_width > unres)
            vel_width = vel_width[ind]
        nlos = np.shape(vel_width)[0]
        #print 'nlos = ',nlos
        v_table = 10**np.arange(0, np.log10(np.max(vel_width)), dv)
        vbin = np.array([(v_table[i]+v_table[i+1])/2. for i in range(0,np.size(v_table)-1)])
        nn = np.histogram(vel_width,v_table)[0] / (1.*nlos)
        width = np.array([v_table[i+1]-v_table[i] for i in range(0,np.size(v_table)-1)])
        vels=nn/width
        return (vbin, vels)

    def get_col_density(self, elem, ion):
        """Get the column density in each pixel for a given species"""
        try:
            rho = self.metals[(elem, ion)][0]
        except KeyError:
            [rho, vel, temp] = self.SPH_Interpolate_metals(elem, ion)
            self.metals[(elem, ion)] = [rho, vel, temp]
        #Size of a bin in physical m
        binsz = self.box/(1.*self.nbins)*self.KPC*(1+self.red)/self.hubble
        #Convert from physical kg/m^2 to atoms/cm^2
        convert = 1./self.PROTONMASS/1e4/self.lines.get_mass(elem)
        return rho*binsz*convert

    def get_vel(self, elem, ion):
        """Get the velocity for a given species, assuming already calculated"""
        return self.metals[(elem, ion)][1]

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
        dX=self.absorption_distance()
        #Total Number of pixels
        rho = self.get_col_density(elem, ion)
        tot_cells = np.size(rho)
        (tot_f_N, NHI_table) = np.histogram(rho,NHI_table)
        tot_f_N=tot_f_N/(width*dX*tot_cells)
        return (center, tot_f_N)

    def absorption_distance(self):
        """Compute X(z), the absorption distance per sightline (eq. 9 of Nagamine et al 2003)
        in dimensionless units."""
        #h * 100 km/s/Mpc in h/s
        h100=3.2407789e-18
        # in cm/s
        light=2.9979e10
        #Units: h/s   s/cm                        kpc/h      cm/kpc
        return h100/light*(1+self.red)**2*self.box*self.KPC

    def rho_crit(self):
        """Get the critical density at z=0 in units of g cm^-3"""
        #H in units of 1/s
        h100=3.2407789e-18*self.hubble
        #G in cm^3 g^-1 s^-2
        grav=6.672e-8
        rho_crit=3*h100**2/(8*math.pi*grav)
        return rho_crit

    def omega_DLA(self, thresh=20.3, elem = "H", ion = 1):
        """Compute Omega_DLA, the sum of the mass in DLAs, divided by the volume of the spectra, divided by the critical density.
            Ω_DLA = m_p * avg. column density / (1+z)^2 / length of column / rho_c
            Note: If we want the neutral gas density rather than the neutral hydrogen density, divide by 0.76,
            the hydrogen mass fraction.
        """
        #Column density of HI in atoms cm^-2 (physical)
        col_den = self.get_col_density(elem, ion)
        if thresh > 0:
            HIden = np.mean(col_den[np.where(col_den > 10**thresh)])
        else:
            HIden = np.mean(col_den)
        #Avg. Column density of HI in g cm^-2 (comoving)
        HIden = self.PROTONMASS * HIden/(1+self.red)**2
        #Length of column in comoving cm
        length = (self.box*self.KPC/self.hubble)
        #Avg density in g/cm^3 (comoving) divided by critical density in g/cm^3
        omega_DLA=HIden/length/self.rho_crit()
        return omega_DLA

    def get_overden(self, thresh = 10**20.3, elem = "H", ion= 1):
        """
        Get an array of spectral pixels which is True where there is a DLA, False otherwise
        """
        col_den = self.get_col_density(elem, ion)
        return np.greater(col_den, thresh)

    def get_spectra_proj_pos(self):
        """Get the position of the spectra projected to their origin"""
        #if np.mean(self.axis) != self.axis[0] or  self.axis[0] != self.axis[-1]:
        #    raise ValueError("Not all spectra are along the same axis")
        axis = self.axis[0]
        if axis == 1:
            spos = self.cofm[:,1:]
        if axis == 2:
            spos = np.vstack([self.cofm[:,0],self.cofm[:,2]]).T
        if axis == 3:
            spos = self.cofm[:,:2]
        return spos

    def _count_modes(self, rbins2, sdist2, pzpos, nspectra):
        """
        Count the autocorrelation modes in each bin.
        This is the slow part of the autocorrelation calculation!
        For two arrays of length L lying at distance 0 from each other
        there will be 2 * 2 * ( L - k )
        pairs at a distance k
        (except for k = 0, when there will be 2L)
        """
        modes = np.zeros(np.size(rbins2)-1)
        pzpos2 = pzpos**2
        lmodes = 4*(self.nbins-np.arange(self.nbins))
        #Zero mode
        lmodes[0]/=2
        #Add modes corresponding to spectra correlated with themselves
        inbins = np.digitize(pzpos2, rbins2)-1
        modes[inbins]+=nspectra*lmodes
        #Count modes: this is the slow part!
        #Squared distance between spectra
        # Count pixels:
        # For two spectra at distance s there will be
        # 4(L-k) pixels at distance sqrt(s^2 + k^2)
        #This call to digitize is still the slow part!
        #One could write an optimised C version of the next three lines
        lmodes *= 2
        for xx in xrange(nspectra):
            for yy in xrange(xx):
                inbins = np.digitize(sdist2[xx,yy]+pzpos2, rbins2)-1
                modes[inbins]+=lmodes
        return modes

    def autocorr(self, thresh=10**20.3, elem = "H", ion = 1, bins=20):
        """
        Compute the autocorrelation function of DLAs along these spectra in logarithmic bins
        Arguments:
            thresh - DLA threshold column density
            elem - element to use
            ion - ionic species to consider
            bins - Number of bins in resulting autocorrelation
        """
        #Autocorrelation bins in log space
        rbins = np.logspace(0,np.log10(self.box * np.sqrt(3)), num=bins+1)
        rbins[0] = 0
        rbins2 = rbins**2
        spos = self.get_spectra_proj_pos()
        sdist2 = (np.square(np.subtract.outer(spos[:,0],spos[:,0]))+np.square(np.subtract.outer(spos[:,1],spos[:,1])))
        #Position of pixel along spectra axis
        pzpos = np.arange(0,self.nbins)*self.box/self.nbins
        rcent = np.array([(rbins[i]+rbins[i+1])/2. for i in xrange(0,bins)])
        auto = np.zeros(bins)
        dlas = self.get_overden(thresh, elem, ion)
        nspectra = np.shape(dlas)[0]
        #Count modes
        modes = self._count_modes(rbins2, sdist2, pzpos,nspectra)
        #Mean value
        nbar = np.mean(dlas)
        dla_ind = np.where(dlas)
        #Array of spectral distance pairs
        #Squared z-distance between two pixels
        pdist2 = np.square(np.subtract.outer(pzpos, pzpos))
        #There is an obvious symmetry, so only do the upper half of the triangle
        for xx in xrange(nspectra):
            dla_x = np.where(dla_ind[0]==xx)
            if np.size(dla_x) == 0:
                continue
            dla_xp = dla_ind[1][dla_x]
            #Do cross-correlation of the same spectra
            #the slightly odd indexing is because numpy complains without it
            dist = pdist2[dla_xp,:][:, dla_xp]
            #And the correlation
            auto += np.histogram(dist, rbins2)[0]
            for yy in xrange(xx):
                #Continue if no DLAs here
                #Pixel pairs where there are DLAs
                dla_y = np.where(dla_ind[0]==yy)
                if np.size(dla_y) == 0:
                    continue
                dla_yp = dla_ind[1][dla_y]
                # Distance between pixels:
                dist2 = pdist2[dla_xp,:][:, dla_yp]+sdist2[xx,yy]
                #And the correlation
                auto += 2*np.histogram(dist2, rbins2)[0]
        ind = np.where(modes == 0)
        modes[ind]=1
        auto=(auto/modes)/nbar**2 - 1
        modes[ind]=0
        return (rcent, modes, auto)
