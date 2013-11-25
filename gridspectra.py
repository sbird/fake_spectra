# -*- coding: utf-8 -*-
"""Class to generate spectra in the positions where there is a DLA, as known from the grid generation."""

import numpy as np
import hdfsim
import h5py
import spectra
import os.path as path

class GridSpectra(spectra.Spectra):
    """Generate metal line spectra from simulation snapshot"""
    def __init__(self,num, base, numlos=5000, res = 1., cdir = None, dla=True, savefile="rand_spectra_DLA.hdf5", savedir=None, gridfile="boxhi_grid_H2.hdf5"):
        #Load halos to push lines through them
        f = hdfsim.get_file(num, base, 0)
        self.box = f["Header"].attrs["BoxSize"]
        f.close()
        if savedir == None:
            savedir = path.join(base,"snapdir_"+str(num).rjust(3,'0'))
        gridfile = path.join(savedir,gridfile)

        self.NumLos = numlos
        #All through y axis
        axis = np.ones(self.NumLos)
        #Load grid positions
        self.dlaind = self._load_dla_index(gridfile, dla)
        #Re-seed for repeatability
        np.random.seed(23)
        cofm = self.get_cofm()
        spectra.Spectra.__init__(self,num, base, cofm, axis, res, cdir, savefile=savefile,savedir=savedir)

        if dla:
            self.replace_not_DLA(10**20.3)
        else:
            self.replace_not_DLA(10**17)
        print "Found DLAs"


    def get_cofm(self, num = None):
        """Find a bunch of sightline positions known to be where a DLA or an LLS is."""
        if num == None:
            num = self.NumLos

        #Get some random indices into the box.
        index = np.random.random_integers(0,np.size(self.dlaind[:,0])-1,num)
        cofm = np.array([self.dlaind[index,0],self.dlaind[index,0],self.dlaind[index,1]]).T
        #Randomize positions within a cell
        cofm[:,1] += self.celsz*(np.random.random_sample(num)-0.5)
        cofm[:,2] += self.celsz*(np.random.random_sample(num)-0.5)
        #Some sightlines could end up being through the same cell, in rare cases.
        #This is only a problem if you want to compare to a quasar survey with pixels large
        #compared to the grid size.
        return cofm

    def _load_dla_index(self, gridfile, dla=True):
        """Load the positions of DLAs or LLS from savefile"""
        #Load the DLA/LLS positions
        f=h5py.File(gridfile,'r')
        grid_file=f["HaloData"]
        ngrid = np.array(grid_file["ngrid"])
        self.celsz = 1.*self.box/ngrid[0]
        grp = f["abslists"]
        #This is needed to make the dimensions right
        if dla:
            ind = (grp["DLA"][0,:],grp["DLA"][1,:],grp["DLA"][2,:])
        else:
            ind = (grp["LLS"][0,:],grp["LLS"][1,:],grp["LLS"][2,:])
        f.close()
        yslab = (ind[1]+0.5)*self.celsz
        zslab = (ind[2]+0.5)*self.celsz
        return np.array((yslab, zslab))

