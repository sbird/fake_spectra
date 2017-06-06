# -*- coding: utf-8 -*-
"""Class to create spectra spaced on a regular grid through the box"""

import numpy as np

from . import abstractsnapshot as absn
from . import spectra

class GriddedSpectra(spectra.Spectra):
    """Generate metal line spectra from simulation snapshot. Default parameters are BOSS DR9"""
    def __init__(self,num, base, nspec=200, res = 90., savefile="gridded_spectra.hdf5", reload_file=True, **kwargs):
        # get box size from file (either HDF5 or BigFile)
        f = absn.AbstractSnapshotFactory(num, base)
        self.box = f.get_header_attr("BoxSize")
        del f
        self.NumLos = nspec*nspec
        #All through y axis
        axis = np.ones(self.NumLos)
        # get position of skewers (on a regular grid)
        cofm = self.get_cofm()
        spectra.Spectra.__init__(self,num, base, cofm=cofm, axis=axis, res=res, savefile=savefile, reload_file=reload_file, **kwargs)

    def get_cofm(self, num = None):
        """Find a bunch more sightlines: should be overriden by child classes"""
        if num is None:
            num = int(np.sqrt(self.NumLos))
        cofm = np.empty([num*num, 3])
        for nn in range(num):
            for mm in range(num):
                cofm[nn*num+mm] = self.box*np.array([nn, nn, mm])/(1.*num)
        return cofm
