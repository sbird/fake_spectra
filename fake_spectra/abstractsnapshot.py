"""Class to abstract the reading of multiple file types for the spectral generator.
Currently supports hdf5 and bigfile."""

import os
import glob
import numpy as np
import h5py

from . import unitsystem

try:
    import bigfile
except ImportError:
    bigfile = False

def AbstractSnapshotFactory(num, base):
    """Function to get a snapshot in whichever format is present"""
    #First try to open it as an HDF5 snapshot
    try:
        return HDF5Snapshot(num, base)
    except IOError:
        if bigfile is False:
            raise IOError("Not an HDF5 snapshot: ", base)
        try:
            return BigFileSnapshot(num, base)
        except (IOError,bigfile.BigFileError):
            raise IOError("Not a bigfile or HDF5 snapshot: ",base)

class AbstractSnapshot(object):
    """A class to abstract a simulation snapshot, so we can use a uniform interface for both HDF5 and bigfile."""
    def __init__(self):
        #Map of (only changed) block names between HDF5 and bigfile snapshots.
        self.hdf_to_bigfile_map = { "Coordinates" : "Position",
                                    "Velocities": "Velocity", "Masses": "Mass",
                                    "NeutralHydrogenAbundance": "NeutralHydrogenFraction",
                                    "GFM_Metallicity": "Metallicity"
                                  }
        #This requires python 2.7
        self.bigfile_to_hdf_map = {v : k for (k,v) in self.hdf_to_bigfile_map.items()}

    def __del__(self):
        try:
            #Note that we will segfault if bigfile is used after close.
            self._f_handle.close()
        except AttributeError:
            pass

    def get_header_attr(self, attr):
        """Return an attribute of the simulation header"""
        attr = self._f_handle["Header"].attrs[attr]
        if isinstance(attr, np.ndarray) and np.size(attr) == 1:
            try:
                return attr[0]
            except IndexError:
                return attr
        return attr

    def get_kernel(self):
        """Get the integer corresponding to a density kernel for each particle.
           The types understood are defined in absorption.h and are currently:
               0 - Top hat kernel
               1 - SPH cubic spline kernel.
            Further additions should match the types defined in allvars.h of MP-Gadget, which are:
                DENSITY_KERNEL_CUBIC_SPLINE = 1,
                DENSITY_KERNEL_QUINTIC_SPLINE = 2,
                DENSITY_KERNEL_QUARTIC_SPLINE = 4,
        """
        return 1

    def get_n_segments(self):
        """Return the number of segments. Number of files on HDF5,
           but may be whatever in convenient for bigfile."""
        raise NotImplementedError

    def get_blocklen(self, part_type, blockname, segment):
        """Get the length of a block"""
        raise NotImplementedError

    def get_data(self, part_type, blockname, segment):
        """Get the data for a particular block, specified
        using either the BigFile names or the HDF5 names.
        Segment: which data segment to load."""
        raise NotImplementedError

    def get_npart(self):
        """Get number of particles in snapshot. Special for HDF5."""
        return self.get_header_attr("TotNumPart")

    def get_omega_baryon(self):
        """Get omega baryon. Special for HDF5."""
        return self.get_header_attr("OmegaBaryon")

    def get_smooth_length(self, part_type, segment):
        """Gets the smoothing length, adjusting the kernel definition to the one we use.
        Special for HDF5, as Arepo stores something else in the SmoothingLength array.
        The kernel is defined so that the smoothing length is 2*h.
        """
        #There is a different kernel definition, as in gadget the kernel goes from 0 to 2,
        #whereas I put it between zero and 1.
        return self.get_data(part_type, "SmoothingLength",segment=segment)/2

    def get_units(self):
        """Get the base scales for the unit system."""
        length = self.get_header_attr("UnitLength_in_cm")
        mass = self.get_header_attr("UnitMass_in_g")
        vel = self.get_header_attr("UnitVelocity_in_cm_per_s")
        units = unitsystem.UnitSystem(UnitLength_in_cm=length, UnitMass_in_g=mass, UnitVelocity_in_cm_per_s=vel)
        return units

    def get_temp(self,part_type, segment,hy_mass=0.76, units=None):
        """Compute temperature (in K) from internal energy.
           Uses: internal energy
                 electron abundance
                 hydrogen mass fraction (0.76)
           Factor to convert U (J/kg) to T (K) : U = N k T / (γ - 1)
           T = U (γ-1) μ m_P / k_B
           where k_B is the Boltzmann constant
           γ is 5/3, the perfect gas constant
           m_P is the proton mass

           μ = 1 / (mean no. molecules per unit atomic weight)
             = 1 / (X + Y /4 + E)
             where E = Ne * X, and Y = (1-X).
             Can neglect metals as they are heavy.
             Leading contribution is from electrons, which is already included
             [+ Z / (12->16)] from metal species
             [+ Z/16*4 ] for OIV from electrons."""
        #convert U (J/kg) to T (K) : U = N k T / (γ - 1)
        #T = U (γ-1) μ m_P / k_B
        #where k_B is the Boltzmann constant
        #γ is 5/3, the perfect gas constant
        #m_P is the proton mass
        #μ is 1 / (mean no. molecules per unit atomic weight) calculated in loop.
        #Internal energy units are 10^-10 erg/g
        if units is None:
            units = self.get_units()
        ienergy = self.get_data(part_type, "InternalEnergy", segment=segment)*units.UnitInternalEnergy_in_cgs
        #Calculate temperature from internal energy and electron abundance
        nelec = self.get_data(part_type, "ElectronAbundance", segment=segment)
        muienergy = 4 / (hy_mass * (3 + 4*nelec) + 1)*ienergy
        #So for T in K, boltzmann in erg/K, internal energy has units of erg/g
        temp = (units.gamma-1) * units.protonmass / units.boltzmann * muienergy
        return temp

class HDF5Snapshot(AbstractSnapshot):
    """Specialised class for loading HDF5 snapshots"""
    def __init__(self, num, base):
        self._files = sorted(self._get_all_files(num, base))
        self._files.reverse()
        self._f_handle = h5py.File(self._files[0], 'r')
        self._handle_num = 0
        AbstractSnapshot.__init__(self)

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
            fnames = glob.glob(os.path.join(fname, "snapshot_"+snap+"*hdf5"))
        if len(fnames) == 0:
            raise IOError("No files found")
        fnames.sort()
        return [fff for fff in fnames if h5py.is_hdf5(fff) ]

    def get_data(self, part_type, blockname, segment):
        """Get the data for a particular particle type.
           Segment: which file to load from."""
        if blockname in self.bigfile_to_hdf_map.keys():
            blockname = self.bigfile_to_hdf_map[blockname]
        if segment < 0:
            def _getone(ff):
                """Get data from one file"""
                fhandle = h5py.File(ff,'r')
                return np.array(fhandle["PartType"+str(part_type)][blockname])
            return np.concatenate([_getone(ff) for ff in self._files])
        if self._handle_num != segment:
            self._f_handle.close()
            self._f_handle = h5py.File(self._files[segment],'r')
            self._handle_num = segment
        return np.array(self._f_handle["PartType"+str(part_type)][blockname])

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

    def get_blocklen(self, part_type, blockname, segment):
        """Get the length of a block"""
        if blockname in self.bigfile_to_hdf_map.keys():
            blockname = self.bigfile_to_hdf_map[blockname]
        if self._handle_num != segment:
            self._f_handle.close()
            self._f_handle = h5py.File(self._files[segment],'r')
            self._handle_num = segment
        return self._f_handle["PartType"+str(part_type)][blockname].len()

    def get_smooth_length(self, part_type, segment):
        """Figures out if the particles are from AREPO or GADGET
        and computes the smoothing length.
        Note the Volume array in HDF5 is comoving and this returns a comoving smoothing length
        The SPH kernel definition used in Gadget (Price 2011: arxiv 1012.1885)
        gives a normalisation so that rho_p = m_p / h^3
        So the smoothing length for Arepo is Volume^{1/3}
        For gadget the kernel is defined so that the smoothing length is 2*h.
        Arguments:
            Baryon particles from a simulation
        Returns:
            Array of smoothing lengths in code units.
        """
        #Are we arepo? If we are a modern version we should have this array.
        try:
            radius = np.power(self.get_data(part_type, "Volume",segment=segment), 1./3)
        except KeyError:
            #If we don't have a Volume array we are gadget, and
            #the SmoothingLength array is actually the smoothing length.
            #There is a different kernel definition, as in gadget the kernel goes from 0 to 2,
            #whereas I put it between zero and 1.
            try:
                radius=self.get_data(part_type, "SmoothingLength",segment=segment)/2
            except KeyError:
                #Very new Arepo has Density and Masses but no volume.
                density = self.get_data(part_type, "Density",segment=segment)
                mass = self.get_data(part_type, "Masses",segment=segment)
                volume = mass/density
                radius = np.power(volume, 1./3)
        return radius

    def get_kernel(self):
        """Get the integer corresponding to a density kernel for each particle.
           The types understood are defined in absorption.h and are currently:
               0 - Top hat kernel (for Arepo)
               1 - SPH cubic spline kernel (for Gadget).
        """
        #We are Arepo.
        if "Volume" in self._f_handle["PartType0"].keys():
            return 0
        else:
            return 1

class BigFileSnapshot(AbstractSnapshot):
    """Specialised class for loading HDF5 snapshots"""
    def __init__(self, num, base):
        fname = base
        snap=str(num).rjust(3,'0')
        new_fname = os.path.join(base, "PART_"+snap)
        #Check for snapshot directory
        if os.path.exists(new_fname):
            fname = new_fname
        self._f_handle = bigfile.BigFile(fname, 'r')
        if "Header" not in self._f_handle.blocks:
            raise IOError("No BigFile snapshot at",new_fname)
        AbstractSnapshot.__init__(self)

    def get_data(self, part_type, blockname, segment):
        """Get the data for a particular block, specified
        using either the BigFile names or the HDF5 names.
        Segment: which data segment to load."""
        if blockname in self.hdf_to_bigfile_map.keys():
            blockname = self.hdf_to_bigfile_map[blockname]
        try:
            (start, end) = self._segment_to_partlist(part_type = part_type, segment=segment)
            return self._f_handle[str(part_type)+"/"+blockname][start:end]
        except bigfile.BigFileError:
            raise KeyError("Not found:"+str(part_type)+"/"+blockname)

    def get_n_segments(self):
        """Return the number of segments. Number of files on HDF5,
           but may be whatever in convenient for bigfile."""
        return int(np.max([1,np.sum(self.get_npart())/(2*256.**3)]))

    def get_blocklen(self, part_type, blockname, segment):
        """Get the length of a block"""
        if blockname in self.hdf_to_bigfile_map.keys():
            blockname = self.hdf_to_bigfile_map[blockname]
        try:
            (start, end) = self._segment_to_partlist(part_type=part_type, segment=segment)
            if end is not None:
                return end - start
            else:
                return self._f_handle[str(part_type)+"/"+blockname].size - start
        except bigfile.BigFileError:
            raise KeyError("Not found:"+str(part_type)+"/"+blockname)

    def _segment_to_partlist(self, part_type, segment):
        """Get the first and last particle in a segment."""
        if segment < 0:
            return (0, None)
        n_segments = self.get_n_segments()
        one_segment = int(self.get_npart()[part_type]//n_segments)
        if segment < n_segments -1:
            return (one_segment*segment, one_segment*(segment+1))
        else:
            return (one_segment*segment, None)

    def get_kernel(self):
        """Get the integer corresponding to a density kernel for each particle.
           Currently only the SPH cubic spline is understood.
            Further additions should match the types defined in allvars.h of MP-Gadget, which are:
                DENSITY_KERNEL_CUBIC_SPLINE = 1,
                DENSITY_KERNEL_QUINTIC_SPLINE = 2,
                DENSITY_KERNEL_QUARTIC_SPLINE = 4,
        """
        try:
            kernel = self.get_header_attr("DensityKernel")
        except KeyError:
            kernel = 1
        #Other types are not yet supported.
        assert kernel == 1
        return kernel
