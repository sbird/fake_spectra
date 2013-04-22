# -*- coding: utf-8 -*-
"""Class to gather and analyse various metal line statistics"""

import numpy as np
import hdfsim
import h5py
import math
import halocat
import os.path as path
import spectra

class HaloSpectra(spectra.Spectra):
    """Generate metal line spectra from simulation snapshot"""
    def __init__(self,num, base, repeat = 3, minpart = 400, res = 1., savefile="halo_spectra.hdf5"):
        self.savefile = path.join(base,"snapdir_"+str(num).rjust(3,'0'),savefile)
        try:
            self.load_savefile(self.savefile)
            #In this case they will have been loaded from the savefile
            cofm = None
            axis = None
        except (IOError, KeyError):
            #Load halos to push lines through them
            f = hdfsim.get_file(num, base, 0)
            self.OmegaM = f["Header"].attrs["Omega0"]
            self.box = f["Header"].attrs["BoxSize"]
            self.npart=f["Header"].attrs["NumPart_Total"]+2**32*f["Header"].attrs["NumPart_Total_HighWord"]
            min_mass = self.min_halo_mass(minpart)
            f.close()
            (ind, self.sub_mass, cofm, self.sub_radii) = halocat.find_wanted_halos(num, base, min_mass)
            self.sub_cofm = cofm
            self.NumLos = np.size(self.sub_mass)
            #All through y axis
            axis = np.ones(self.NumLos)
            #axis[self.NumLos/3:2*self.NumLos/3] = 2
            #axis[2*self.NumLos/3:self.NumLos] = 3
            axis = np.repeat(axis,repeat)
            self.NumLos*=repeat
            self.repeat = repeat
            #Re-seed for repeatability
            np.random.seed(23)
            cofm = self.get_cofm()

        spectra.Spectra.__init__(self,num, base, cofm, axis, res, savefile=self.savefile)

        self.replace_not_DLA()

    def replace_not_DLA(self):
        """
        Replace those sightlines which do not contain a DLA with new sightlines, until all sightlines contain a DLA.
        """
        ind = self.filter_DLA()
        while np.size(ind) > 0:
            #Replace spectra that did not result in a DLA
            cofm_new = self.get_cofm()
            self.cofm[ind] = cofm_new[ind]
            [rho, vel, temp] = self.SPH_Interpolate_metals("H", 1, ind = ind)
            self.metals[("H", 1)][0][ind] = rho
            self.metals[("H", 1)][1][ind] = vel
            self.metals[("H", 1)][2][ind] = temp
            ind = self.filter_DLA()

    def get_cofm(self):
        """Find a bunch more sightlines"""
        cofm = np.repeat(self.sub_cofm,self.repeat,axis=0)
        #Perturb the sightlines within a sphere of the virial radius.
        #We want a representative sample of DLAs.
        maxr = self.sub_radii
        #Generate random sphericals
        theta = 2*math.pi*np.random.random_sample(self.NumLos)-math.pi
        phi = 2*math.pi*np.random.random_sample(self.NumLos)
        rr = np.repeat(maxr,self.repeat)*np.random.random_sample(self.NumLos)
        #Add them to halo centers
        cofm[:,0]+=rr*np.sin(theta)*np.cos(phi)
        cofm[:,1]+=rr*np.sin(theta)*np.sin(phi)
        cofm[:,2]+=rr*np.cos(theta)
        return cofm

    def filter_DLA(self):
        """Find sightlines without a DLA"""
        col_den = self.get_col_density("H",1)
        ind = np.where(np.max(col_den, axis=1) < 10**20.3)
        return ind

    def save_file(self):
        """
        Save additional halo data to the savefile
        """
        try:
            f=h5py.File(self.savefile,'a')
        except IOError:
            raise IOError("Could not open ",self.savefile," for writing")
        grp = f.create_group("halos")
        grp["radii"] = self.sub_radii
        grp["cofm"] = self.sub_cofm
        grp["mass"] = self.sub_mass
        grp.attrs["repeat"] = self.repeat
        grp.attrs["NumLos"] = self.NumLos
        f.close()
        spectra.Spectra.save_file(self)

    def load_savefile(self,savefile=None):
        """Load data from a file"""
        #Name of savefile
        f=h5py.File(savefile,'r')
        grp = f["halos"]
        self.sub_radii = np.array(grp["radii"])
        self.sub_cofm = np.array(grp["cofm"])
        self.sub_mass = np.array(grp["mass"])
        self.repeat = grp.attrs["repeat"]
        self.NumLos = grp.attrs["NumLos"]
        f.close()
        spectra.Spectra.load_savefile(self, savefile)


    def min_halo_mass(self, minpart = 400):
        """Min resolved halo mass in internal Gadget units (1e10 M_sun)"""
        #This is rho_c in units of h^-1 1e10 M_sun (kpc/h)^-3
        rhom = 2.78e+11* self.OmegaM / 1e10 / (1e3**3)
        #Mass of an SPH particle, in units of 1e10 M_sun, x omega_m/ omega_b.
        target_mass = self.box**3 * rhom / self.npart[0]
        min_mass = target_mass * minpart
        return min_mass

    def find_associated_halo(self, num):
        """Find the halo sightline num is associated with"""
        nh = num /self.repeat
        return (nh, self.sub_mass[nh], self.sub_cofm[nh,:], self.sub_radii[nh])

    def line_offsets(self):
        """Find the distance between each line and its parent halo"""
        offsets = np.zeros(self.NumLos)
        for ii in xrange(self.NumLos):
            offsets[ii] = np.sqrt(np.sum((self.find_associated_halo(ii)[2]-self.cofm[ii])**2))
        return offsets
