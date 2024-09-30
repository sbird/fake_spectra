# -*- coding: utf-8 -*-
"""Class to create spectra spaced on a regular grid through the box"""

import numpy as np

from . import abstractsnapshot as absn
from . import spectra

class GriddedSpectra(spectra.Spectra):
    """Generate regular grid of spectra along a given axis."""

    def __init__(self,num, base, nspec=200, MPI=None, res = None, nbins=None, savefile="gridded_spectra.hdf5", reload_file=True, grid_axes=None, grid_cofm=None, axis=1, **kwargs):
        """
        Parameters
        ----------
        num : int, required
            Snapshot number
        base : str, required
            Base directory of the snapshot, i.e. containing .hdf5 or PART_XXXX files
        nspec : int, optional
            Number of skewers per side of the grid in the transverse direction
        MPI : object, optional
            MPI communicator
        res : float, required if `nbins` is not provided
            Resolution of the skewers in km/s
        nbins : int, required if `res` is not provided
            Number of bins in the skewers
        savefile : str, optional
            Name of the file to save the spectra
        reload_file : bool, optional
            Try to load the spectra file if it exists
        grid_axes : array of shape (nspec*nspec*nspec,3), optional
            Axes of the grid. Onlly if not a unifrom grid across the box is desired.
            if none, `get_axes_and_cofm` will be called to get the axes
        grid_cofm : array of shape (nspec*nspec*nspec,3), optional
            Centers of the grid. Onlly if not a unifrom grid across the box is desired
            if none, `get_axes_and_cofm` will be called to get the centers
        axis : int, optional
            Axis along which to place the skewers. 
            Only if `grid_axes` and `grid_cofm` are not provided.
        """

        if reload_file is False:
            # if reading skewers from file, no need to read the snapshot
            grid_cofm=None
            grid_axes=None
        else:
            # get box size from file (either HDF5 or BigFile)
            f = absn.AbstractSnapshotFactory(num, base, MPI=MPI)
            self.box = f.get_header_attr("BoxSize")
            del f
            # get position of skewers in the grid
            if grid_cofm is None:
                grid_axes, grid_cofm = self.get_axes_and_cofm(nspec,axis)

        # call constructor of base class
        spectra.Spectra.__init__(self,num,base,cofm=grid_cofm,axis=grid_axes,MPI=MPI, nbins=nbins, res=res,savefile=savefile,reload_file=reload_file,**kwargs)


    def get_axes_and_cofm(self,nspec,axis):
        """Define position of skewers in the grid"""

        # decide axis to use to extract skewers (x axis by default)
        if axis < 0:
            grid_axes = np.ones(3*nspec*nspec)
            grid_axes[nspec*nspec:2*nspec*nspec] = 2
            grid_axes[2*nspec*nspec:] = 3
        else:
            grid_axes = axis*np.ones(nspec*nspec)

        # figure out grid of positions for this particular axis
        if axis==1:
            grid_id=np.array([np.array([0,nn,mm]) for nn in range(nspec) for mm in range(nspec)])
        elif axis==2:
            grid_id=np.array([np.array([nn,0,mm]) for nn in range(nspec) for mm in range(nspec)])
        elif axis==3:
            grid_id=np.array([np.array([nn,mm,0]) for nn in range(nspec) for mm in range(nspec)])
        elif axis < 0:
            grid_id_1=np.array([np.array([0,nn,mm]) for nn in range(nspec) for mm in range(nspec)])
            grid_id_2=np.array([np.array([nn,0,mm]) for nn in range(nspec) for mm in range(nspec)])
            grid_id_3=np.array([np.array([nn,mm,0]) for nn in range(nspec) for mm in range(nspec)])
            grid_id=np.concatenate([grid_id_1, grid_id_2, grid_id_3])
        else:
            raise ValueError('wrong axis number {}'.format(axis))

        # separation between skewers in the grid (per side)
        dx=self.box/(1.*nspec)

        return grid_axes, dx*grid_id
