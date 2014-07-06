# -*- coding: utf-8 -*-
"""Class to generate spectra in the positions where there is a DLA, as known from the grid generation."""

import numpy as np
import hdfsim
import h5py
import vw_spectra
import os.path as path

class GridSpectra(vw_spectra.VWSpectra):
    """Generate metal line spectra from simulation snapshot"""
    def __init__(self,num, base, numlos=5000, res = 1., cdir = None, dla=True, savefile="grid_spectra_DLA.hdf5", savedir=None, gridfile="boxhi_grid_H2.hdf5"):
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
        self.dlaind = self._load_dla_index(gridfile)
        #Re-seed for repeatability
        np.random.seed(23)
        cofm = self.get_cofm()
        spectra.Spectra.__init__(self,num, base, cofm, axis, res, cdir, savefile=savefile,savedir=savedir, reload_file=True)

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
        index = np.random.random_integers(0,np.size(self.dlaind[0,:])-1,num)
        cofm = np.array([self.dlaind[0,index],self.dlaind[0,index],self.dlaind[1,index]]).T
        #Randomize positions within a cell
        cofm[:,1] += self.celsz*(np.random.random_sample(num)-0.5)
        cofm[:,2] += self.celsz*(np.random.random_sample(num)-0.5)
        #Some sightlines could end up being through the same cell, in rare cases.
        #This is only a problem if you want to compare to a quasar survey with pixels large
        #compared to the grid size.
        return cofm

    def _load_dla_index(self, gridfile, dla=False):
        """Load the positions of DLAs or LLS from savefile"""
        #Load the DLA/LLS positions
        f=h5py.File(gridfile,'r')
        grid_file=f["HaloData"]
        ngrid = np.array(grid_file["ngrid"])
        self.celsz = 1.*self.box/ngrid[0]
        grp = f["abslists"]
        #This is needed to make the dimensions right
        ind = (grp["DLA"][0,:],grp["DLA"][1,:],grp["DLA"][2,:])
        if not dla:
            ind_lls = (grp["LLS"][0,:],grp["LLS"][1,:],grp["LLS"][2,:])
        f.close()
        yslab = (ind[1]+0.5)*self.celsz
        yslab_lls = (ind_lls[1]+0.5)*self.celsz
        yslab = np.append(yslab,yslab_lls)
        zslab = (ind[2]+0.5)*self.celsz
        zslab_lls = (ind_lls[2]+0.5)*self.celsz
        zslab = np.append(zslab,zslab_lls)
        return np.array((yslab, zslab))


class TestGridSpectra(GridSpectra):
    """This specialised class tests the spectral generation code by loading several sightlines in a single cell and finding
    their average value, compared to the value in the cell."""
    def __init__(self,num, base, numlos=5000, res = 1., seed=23,cdir = None, dla=True, savefile="grid_spectra_DLA.hdf5", savedir=None, gridfile="boxhi_grid_H2.hdf5"):
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
        self.dlaval = self._load_dla_val(gridfile, dla)
        #Re-seed for repeatability
        np.random.seed(seed)
        cofm = self.get_cofm()
        spectra.Spectra.__init__(self,num, base, cofm, axis, res, cdir, savefile=savefile,savedir=savedir, reload_file=True)

    def get_cofm(self, num = None):
        """Find a bunch of sightline positions through a single cell containing a DLA."""
        if num == None:
            num = self.NumLos

        #Get a single random position
        self.index = np.random.random_integers(0,np.size(self.dlaind[0,:])-1,1)*np.ones(num,dtype=np.int)
        cofm = np.array([self.dlaind[0,self.index],self.dlaind[0,self.index],self.dlaind[1,self.index]]).T
        #Randomize positions within a cell
        cofm[:,1] += self.celsz*(np.random.random_sample(num)-0.5)
        cofm[:,2] += self.celsz*(np.random.random_sample(num)-0.5)
        #Some sightlines could end up being through the same cell, in rare cases.
        #This is only a problem if you want to compare to a quasar survey with pixels large
        #compared to the grid size.
        return cofm

    def check_mean(self):
        """Compute difference between the mean column of the spectra in this cell and the grid value."""
        dlaval = self.dlaval[self.index][0]
        colden = self.get_col_density("H",1)
        specval = np.sum(colden)/self.NumLos
        print "From spectra:",specval
        print "From grid:",10**dlaval
        print "different:",specval/10**dlaval

    def _load_dla_val(self, gridfile, dla=True):
        """Load the values of DLAs or LLS from savefile"""
        #Load the DLA/LLS positions
        f=h5py.File(gridfile,'r')
        grp = f["abslists"]
        #This is needed to make the dimensions right
        if dla:
            nhi = np.array(grp["DLA_val"])
        else:
            nhi = np.array(grp["LLS_val"])
        f.close()
        return nhi
