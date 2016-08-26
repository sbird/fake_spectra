# -*- coding: utf-8 -*-
"""Class to create spectra spaced on a regular grid through the box"""

import numpy as np
import hdfsim
import spectra

class GriddedSpectra(spectra.Spectra):
    """Generate metal line spectra from simulation snapshot. Default parameters are BOSS DR9"""
    def __init__(self,num, base, nspec=200, res = 90., cdir = None, savefile="gridded_spectra.hdf5", savedir=None):
        #Load halos to push lines through them
        f = hdfsim.get_file(num, base, 0)
        self.box = f["Header"].attrs["BoxSize"]
        f.close()
        self.NumLos = nspec*nspec
        #All through y axis
        axis = np.ones(self.NumLos)
        #Sightlines at random positions
        #Re-seed for repeatability
        np.random.seed(23)
        cofm = self.get_cofm()
        spectra.Spectra.__init__(self,num, base, cofm, axis, res, cdir, savefile=savefile,savedir=savedir,reload_file=True)

    def get_cofm(self, num = None):
        """Find a bunch more sightlines: should be overriden by child classes"""
        if num is None:
            num = int(np.sqrt(self.NumLos))
        cofm = np.empty([num*num, 3])
        for nn in range(num):
            for mm in range(num):
                cofm[nn*num+mm] = self.box*np.array([nn, nn, mm])/(1.*num)
        return cofm
