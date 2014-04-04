# -*- coding: utf-8 -*-
"""Class to gather and analyse various metal line statistics"""

import numpy as np
import hdfsim
import spectra

class CIVSpectra(spectra.Spectra):
    """Generate metal line spectra from simulation snapshot"""
    def __init__(self,num, base, numlos=5000, res = 1., cdir = None, thresh=0.05, savefile="civ_spectra_DLA.hdf5", savedir=None):
        #Load halos to push lines through them
        f = hdfsim.get_file(num, base, 0)
        self.box = f["Header"].attrs["BoxSize"]
        f.close()
        self.NumLos = numlos
        #All through y axis
        axis = np.ones(self.NumLos)
        #Sightlines at random positions
        #Re-seed for repeatability
        np.random.seed(23)
        cofm = self.get_cofm()
        spectra.Spectra.__init__(self,num, base, cofm, axis, res, cdir, savefile=savefile,savedir=savedir,reload_file=True)

        if thresh > 0:
            self.replace_not_DLA(thresh)
        print "Found DLAs"


    def get_cofm(self, num = None):
        """Find a bunch more sightlines: should be overriden by child classes"""
        if num == None:
            num = self.NumLos
        cofm = self.box*np.random.random_sample((num,3))
        return cofm

    def replace_not_DLA(self, thresh=0.05):
        """
        Replace those sightlines which do not contain CIV with eq. width > thresh with new sightlines.
        Must implement get_cofm for this to work
        """
        #Declare variables
        found = 0
        wanted = self.NumLos
        cofm_DLA = np.empty_like(self.cofm)
        #Filter
        eqw = self.equivalent_width("C",4,1550)
        del self.tau[("C",4,1550)]
        ind = np.where(eqw > thresh)
        #Update saves
        top = np.min([wanted, found+np.size(ind)])
        cofm_DLA[found:top] = self.cofm[ind][:top,:]
        found += np.size(ind)
        self.discarded = wanted-np.size(ind)
        print "Discarded: ",self.discarded
        while found < wanted:
            #Get a bunch of new spectra
            self.cofm = self.get_cofm()
            eqw = self.equivalent_width("C",4,1550)
            ind = np.where(eqw > thresh)
            del self.tau[("C",4,1550)]
            #Update saves
            top = np.min([wanted, found+np.size(ind)])
            cofm_DLA[found:top] = self.cofm[ind][:top-found,:]
            found += np.size(ind)
            self.discarded += wanted-np.size(ind)
            print "Discarded: ",self.discarded

        #Copy back
        self.cofm=cofm_DLA
        #Finalise the cofm array
        self.cofm_final = True

