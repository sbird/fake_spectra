"""Class to abstract the reading of multiple file types for the spectral generator.
Currently supports hdf5 and bigfile."""

import os
import glob
import numpy as np
import h5py
try:
    import bigfile
except ImportError:
    pass

def AbstractSnapshotFactory(num, base):
    """Function to get a snapshot in whichever format is present"""
    #First try to open it as an HDF5 snapshot
    try:
        return HDF5Snapshot(num, base)
    except IOError:
        raise IOError("Not an HDF5 snapshot set")

class AbstractSnapshot(object):
    """A class to abstract a simulation snapshot, so we can use a uniform interface for both HDF5 and bigfile."""
    _f_handle = None
    def __del__(self):
        self._f_handle.close()
        #Note that we will segfault if bigfile is used after close.

    def get_header_attr(self, attr):
        """Return an attribute of the simulation header"""
        return self._f_handle["Header"].attrs[attr]

    def get_n_segments(self):
        """Return the number of segments. Number of files on HDF5,
           but may be whatever in convenient for bigfile."""
        raise NotImplementedError

    def get_particle_data(self, part_type, segment=0):
        """Get the data for a particular particle type.
        Segment: which data segment to load"""
        raise NotImplementedError

    def get_npart(self):
        """Get number of particles in snapshot. Special for HDF5."""
        return self.get_header_attr("TotNumPart")

    def get_omega_baryon(self):
        """Get omega baryon. Special for HDF5."""
        return self.get_header_attr("OmegaBaryon")

class HDF5Snapshot(AbstractSnapshot):
    """Specialised class for loading HDF5 snapshots"""
    def __init__(self, num, base):
        self._files = sorted(self._get_all_files(num, base))
        self._files.reverse()
        self._f_handle = h5py.File(self._files[0], 'r')
        self._handle_num = 0

    def _get_all_files(self, num, base):
        """Get a file descriptor from a simulation directory,
        snapshot number and optionally file number.
        Input:
            num - snapshot number
            base - simulation directory
            file_num - file number in the snapshot"""
        fname = base
        snap=str(num).rjust(3,'0')
        new_fname = os.path.join(base, "snapdir_"+snap)
        #Check for snapshot directory
        if os.path.exists(new_fname):
            fname = new_fname
        #Find a file
        fnames = glob.glob(os.path.join(fname, "snap_"+snap+"*hdf5"))
        if len(fnames) == 0:
            raise IOError("No files found")
        fnames.sort()
        return [fff for fff in fnames if h5py.is_hdf5(fff) ]

    def get_particle_data(self, part_type, segment=0):
        """Get the data for a particular particle type.
           Segment: which file to load from."""
        if self._handle_num != segment:
            self._f_handle.close()
            self._f_handle = h5py.File(self._files[segment],'r')
            self._handle_num = segment
        return self._f_handle["PartType"+str(part_type)]

    def get_npart(self):
        """Get the total number of particles in the snapshot."""
        return self.get_header_attr("NumPart_Total")+2**32*self.get_header_attr("NumPart_Total_HighWord")

    def get_omega_baryon(self):
        """Get omega_baryon from a single file (neglecting stars)."""
        mass_dm = self.get_header_attr("MassTable")[1]*self.get_header_attr("NumPart_ThisFile")[1]
        mass_bar = np.sum(self._f_handle["PartType0"]["Masses"])
        return mass_bar/(mass_bar+mass_dm)*self.get_header_attr("Omega0")

    def get_n_segments(self):
        """Return the number of segments. Number of files on HDF5,
           but may be whatever in convenient for bigfile."""
        return len(self._files)
