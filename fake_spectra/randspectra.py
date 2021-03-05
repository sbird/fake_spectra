# -*- coding: utf-8 -*-
"""Class to gather and analyse various metal line statistics"""

from __future__ import print_function
import numpy as np

from . import abstractsnapshot as absn
from . import spectra

class RandSpectra(spectra.Spectra):
    """Generate metal line spectra from simulation snapshot"""
    def __init__(self, num, base,  MPI, comm, seed=23,ndla = 1000, numlos=5000, thresh=10**20.3, savefile="rand_spectra_DLA.hdf5", elem="H", ion=1,**kwargs):
        #Load halos to push lines through them
        f = absn.AbstractSnapshotFactory(num, base, comm)
        self.box = f.get_header_attr("BoxSize")
        del f
        self.NumLos = numlos
        #All through z axis (Fortran convention, 1 for x, 2 for y, 3 for z)
        axis = np.ones(self.NumLos)*3
        #Sightlines at random positions
        #Re-seed for repeatability
        np.random.seed(seed)
        cofm = self.get_cofm()
        spectra.Spectra.__init__(self,num, base, MPI, comm, cofm, axis, savefile=savefile,reload_file=True, load_halo=False, **kwargs)
        
        ### I will call replace_not_DLA, in parallel script (After adding col_den of all hdf5 files)
        
        if np.size(thresh) > 1 or thresh > 0:
            #self.initial_spectra(ndla, thresh, elem=elem, ion=ion)
            self.replace_not_DLA(ndla, thresh, elem=elem, ion=ion)
            print("Found objects over threshold")
        
    def get_cofm(self, num = None):
        """Find a bunch more sightlines: should be overriden by child classes"""
        if num is None:
            num = self.NumLos
        cofm = self.box*np.random.random_sample((num,3))
        return cofm
