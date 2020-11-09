# -*- coding: utf-8 -*-
"""Class to gather and analyse various metal line statistics"""

import math
import os.path as path
import numpy as np
import h5py

from . import spectra
from . import halocat

try:
    xrange(1)
except NameError:
    xrange = range

class HaloSpectra(spectra.Spectra):
    """Generate spectra from simulation snapshot which are near to galactic halos."""
    def __init__(self,num, base, repeat = 1, min_mass = 1e9, max_mass=1e11, res = None, offset=1., savefile="halo_spectra.hdf5", savedir=None, cdir=None):
        if savedir is None:
            savedir = path.join(base,"snapdir_"+str(num).rjust(3,'0'))
        self.savefile = path.join(savedir,savefile)
        try:
            self.load_savefile(self.savefile)
            #In this case they will have been loaded from the savefile
            cofm = None
            axis = None
        except (IOError, KeyError):
            #Load halos to push lines through them
            min_mass = min_mass/1e10
            (_, self.sub_mass, cofm, self.sub_radii) = halocat.find_wanted_halos(num, base, min_mass)
            ii = np.where(self.sub_mass < max_mass)
            self.sub_mass = self.sub_mass[ii]
            self.sub_radii = self.sub_radii[ii]
            self.sub_cofm = np.array(cofm[ii], dtype=np.float64)
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
            cofm = self.get_cofm(offset)

        spectra.Spectra.__init__(self,num, base, cofm, axis, res, savefile=self.savefile, savedir=savedir,reload_file=True, cdir=cdir)

        #If we did not load from a snapshot
        if np.size(cofm) > 0:
            self.replace_not_DLA(self.NumLos)

    def get_cofm(self, offset=1, *, num = None):
        """Find a bunch more sightlines"""
        if num != None:
            raise NotImplementedError
        cofm = np.repeat(self.sub_cofm,self.repeat,axis=0)
        #Perturb the sightlines within a sphere of the virial radius.
        maxr = self.sub_radii*offset
        #Generate random sphericals
        theta = 2*math.pi*np.random.random_sample(self.NumLos)-math.pi
        phi = 2*math.pi*np.random.random_sample(self.NumLos)
        rr = np.repeat(maxr,self.repeat)*np.random.random_sample(self.NumLos)
        #Add them to halo centers
        cofm[:,0]+=rr*np.sin(theta)*np.cos(phi)
        cofm[:,1]+=rr*np.sin(theta)*np.sin(phi)
        cofm[:,2]+=rr*np.cos(theta)
        return cofm

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


    def find_associated_halo(self, num):
        """Find the halo sightline num is associated with"""
        nh = num /self.repeat
        return (nh, self.sub_mass[nh], self.sub_cofm[nh,:], self.sub_radii[nh])

    def line_offsets(self):
        """Find the minimum distance between each line and its parent halo"""
        offsets = np.zeros(self.NumLos)
        hcofm = np.array([self.find_associated_halo(ii)[2] for ii in range(self.NumLos)])
        hrad = np.array([self.find_associated_halo(ii)[3] for ii in range(self.NumLos)])
        hpos = self.get_spectra_proj_pos(cofm=hcofm)
        lpos = self.get_spectra_proj_pos()
        offsets = np.sqrt(np.sum((hpos-lpos)**2,axis=1))/hrad
        return offsets

    def load_halo(self):
        """Do nothing - halos already loaded"""
        return

    def replace_not_DLA(self, ndla, thresh=10**20.3, elem="H", ion=1):
        """Do nothing"""
        return
