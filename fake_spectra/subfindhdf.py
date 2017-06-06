"""Small class for reading in subfind tables in HDF5 format"""
import glob
import os.path as path
import h5py
import numpy as np

class SubFindHDF5:
    """
    Class to read the subfind tables from a directory.
    Reads all arrays from all tables and concatenates them together.
    """
    def __init__(self, base, num):
        snap=str(num).rjust(3,'0')
        fofdir = path.join(base,"groups_"+snap)
        fofpatt = fofdir+"/fof_subhalo_tab_"+snap+"*.hdf5"
        #Get a list of the files
        self.foffiles = sorted(glob.glob(fofpatt))
        f = h5py.File(self.foffiles[0],'r')
        #Find out how many of each type of thing we have by reading the header
        self._sizes = {}
        self._sizes["Group"] = f["Header"].attrs["Ngroups_Total"]
        self._sizes["Subhalo"] = f["Header"].attrs["Nsubgroups_Total"]
        #This is not actually used yet
        self._sizes["IDs"] = f["Header"].attrs["Nids_Total"]
        #Find the group and subhalo array names
        self.Grpnames = list(f["Group"].keys())
        self.Subnames = list(f["Subhalo"].keys())
        f.close()
        self._cache = {}
        self._cache["Group"] = {}
        self._cache["Subhalo"] = {}

    def get_grp_names(self):
        """Get the names of all arrays attached to groups"""
        return self.Grpnames

    def get_sub_names(self):
        """Get the names of all arrays attached to subhalos"""
        return self.Subnames

    def _get_single_file_array(self, fname, dset, name, shapes):
        """Get the desired dataset from a single file"""
        f = h5py.File(fname,'r')
        try:
            tmp = np.array(f[dset][name])
        except KeyError:
            #We want an empty array which fits with the full one
            if len(shapes) > 1:
                shapes_zero = (0,)+shapes[1:]
                tmp = np.array([]).reshape(shapes_zero)
            else:
                tmp = np.array([])
        finally:
            f.close()
        return tmp

    def _get_array(self, dset, name):
        """Get the array called 'name' from the array 'dset'"""
        try:
            return self._cache[dset][name]
        except KeyError:
            #Get the shape from the first file, which will always contain a halo.
            shapes = np.shape(self._get_single_file_array(self.foffiles[0], dset, name, None))
            data = np.concatenate([self._get_single_file_array(ff, dset, name, shapes) for ff in self.foffiles])
        #Check we found everything
        assert np.shape(data)[0] == self._sizes[dset]
        self._cache[dset][name] = data
        return data

    def get_sub(self, name):
        """Get the array called 'name' from the Subhalo arrays"""
        return self._get_array("Subhalo", name)

    def get_grp(self, name):
        """Get the array called 'name' from the Group arrays"""
        return self._get_array("Group", name)
