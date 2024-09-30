"""Class to abstract the reading of multiple file types for the spectral generator.
Currently supports hdf5 and bigfile."""

import os
import glob
import numpy as np
import h5py
import logging

from . import unitsystem

def AbstractSnapshotFactory(num, base, MPI=None, log_level='info'):
    """Function to get a snapshot in whichever format is present"""
    #First try to open it as an HDF5 snapshot
    try:
        return HDF5Snapshot(num, base, MPI)
    except IOError:
        try:
            return BigFileSnapshot(num, base, MPI, log_level)
        except:
            raise IOError("Not a bigfile or HDF5 snapshot: ",base)

class AbstractSnapshot(object):
    """A class to abstract a simulation snapshot, so we can use a uniform interface for both HDF5 and bigfile."""
    def __init__(self, log_level='info'):
        #Map of (only changed) block names between HDF5 and bigfile snapshots.
        self.hdf_to_bigfile_map = { "Coordinates" : "Position",
                                    "Velocities": "Velocity", "Masses": "Mass",
                                    "NeutralHydrogenAbundance": "NeutralHydrogenFraction",
                                    "GFM_Metallicity": "Metallicity",
                                    "GFM_Metals": "Metals"
                                  }
        #This requires python 2.7
        self.bigfile_to_hdf_map = {v : k for (k,v) in self.hdf_to_bigfile_map.items()}
        #Set up logging
        logging_level = getattr(logging, log_level.upper(), None)
        self.configure_logging(logging_level)

    def __del__(self):
        try:
            #Note that we will segfault if bigfile is used after close.
            self._f_handle.close()
        except AttributeError:
            pass

    def configure_logging(self, logging_level):
        """Sets up logging based on the provided logging level."""
        self.logger = logging.getLogger('MockAbsorp')
        self.logger.setLevel(logging_level)

        console_handler = logging.StreamHandler()
        formatter = logging.Formatter('%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
        console_handler.setFormatter(formatter)
        self.logger.addHandler(console_handler)
        if self.rank==0:
            self.logger.debug('Logger initialized at level: %s', logging.getLevelName(logging_level))

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
        """Get number of particles in snapshot."""
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
        try:
            length = self.get_header_attr("UnitLength_in_cm")
            mass = self.get_header_attr("UnitMass_in_g")
            vel = self.get_header_attr("UnitVelocity_in_cm_per_s")
            units = unitsystem.UnitSystem(UnitLength_in_cm=length, UnitMass_in_g=mass, UnitVelocity_in_cm_per_s=vel)
        except KeyError:
            print("Warning: Using default (kpc,10^10Msun, km/s) units.")
            units = unitsystem.UnitSystem()
        return units

    def get_peculiar_velocity(self, part_type, segment):
        """Get the peculiar velocity in internal units. Converts out the various Gadget comoving a factors."""
        vel = self.get_data(part_type, "Velocities", segment = segment)
        atime = self.get_header_attr("Time")
        vel *= np.sqrt(atime)
        return vel

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
    def __init__(self, num, base, MPI, log_level='info'):
        #MPI communicator
        self.MPI = MPI
        self.comm = self.MPI.COMM_WORLD
        self._files = sorted(self._get_all_files(num, base))
        self._files.reverse()
        self._f_handle = h5py.File(self._files[0], 'r')
        self._handle_num = 0
        AbstractSnapshot.__init__(self, log_level=log_level)

    def _get_all_files(self, num, base):
        """Get a file descriptor from a simulation directory,
        snapshot number and optionally file number.
        Input:
            num - snapshot number
            base - simulation directory
            file_num - file number in the snapshot"""
        fname = base
        snap=str(num).rjust(3,'0')

        if self.comm is not None:
            rank = self.comm.Get_rank()
            size = self.comm.Get_size()
            new_fname = base
        else :
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
        
        if self.comm is None:
            return [fff for fff in fnames if h5py.is_hdf5(fff) ]
        
        else:
            num_files = len(fnames)
            files_per_rank = int(num_files/size)
            # a list if file names for each rank
            fnames_rank = fnames[rank*files_per_rank : (rank+1)*files_per_rank]
            #some ranks get 1 more snapshot file
            remained = int(num_files - (files_per_rank*size))
            if rank in range(1,remained+1):
                fnames_rank.append(fnames[files_per_rank*size + rank-1 ])
            return [fff for fff in fnames_rank if h5py.is_hdf5(fff) ]

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

    def get_n_segments(self, part_type=None):
        """Return the number of segments. Number of files on HDF5,
           but may be whatever in convenient for bigfile. part_type
           is not used here, it is gotten just for compatibility.
        """
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
               2 - A partial reconstruction of the Voronoi mesh, for Arepo.
               3 - SPH quintic spline kernel (for modern SPH Gadget)
        """
        #We are an older Arepo version if there is a Volume key
        if "Volume" in self._f_handle["PartType0"].keys():
            return 0
        #We are Gadget
        if "SmoothingLength" in self._f_handle["PartType0"].keys():
            return 1
        #We are Arepo if there is no smoothing length.
        #return 0 for the tophat kernel. 2 is also an option,
        #but has too many edge cases to be the default.
        return 0

class BigFileSnapshot(AbstractSnapshot):
    """Specialised class for loading bigfile snapshots"""
    def __init__(self, num, base, MPI=None, log_level='info'):
        self.MPI = MPI
        if self.MPI is not None :
            self.comm = self.MPI.COMM_WORLD
            self.rank = self.comm.Get_rank()
            self.size = self.comm.Get_size()
        else :
            self.comm = None
            self.size = 1
            self.rank = 0
        AbstractSnapshot.__init__(self, log_level=log_level)
        self.parts_rank = None
        # MPI load on each rank to be computed
        self.offset = None
        self.chunk_size = None
        
        part_dir = base
        snap=str(num).rjust(3,'0')
        new_part_dir = os.path.join(part_dir, "PART_"+snap)
        #Check for snapshot directory
        if os.path.exists(new_part_dir):
            part_dir = new_part_dir 
        self.part_dir = part_dir
        # load the snapshot header on root rank and broadcast it
        header = None
        if self.rank == 0:
            header = self._load_header(part_dir)
        if self.comm is not None:
            header = self.comm.bcast(header, root=0)
        self.header = header
        # Placeholfer to store the block headers
        self.block_headers = {}
    
    def _load_header(self, fname):
        """Load the header of the BigFile snapshot
        To be called by only the root rank
        Note: We don't use bigfile API for now to avoid MPI issues
        with the python wrapper.
        Parameters:
            fname (str): The path to the BigFile snapshot
        """
        def process_line(line):
            # Split the line into the main content and the human-readable part
            line_parts = line.split('#HUMANE')
            first_part = line_parts[0].strip()
            second_part = line_parts[1].strip()

            # Split the main content to extract key, data type, size, and hexadecimal data
            first_parts = first_part.split()
            key = first_parts[0]
            dtype_str = first_parts[1]
            size = int(first_parts[2])
            hex_data = first_parts[3]

            # Convert the hexadecimal data to bytes
            data_bytes = bytes.fromhex(hex_data)

            # Process the data based on its type
            if dtype_str.startswith('<S'):
                # For strings, read the data as a single string of the specified size
                dtype = np.dtype('S{}'.format(size))
                data = np.frombuffer(data_bytes, dtype=dtype, count=1)[0].decode('ascii').rstrip('\x00')
            else:
                # For numerical data, read the data according to its type and size
                dtype = np.dtype(dtype_str)
                data = np.frombuffer(data_bytes, dtype=dtype, count=size)
                # If there's only one element, extract it directly
                if size == 1:
                    data = data[0]
                else:
                    data = data.tolist()

            return key, data
        header = {}
        with open(os.path.join(fname, 'Header','attr-v2')) as f:
            for line in f:
                line = line.strip()
                if line:
                    key, data = process_line(line)
                    header[key] = data
        return header         

    def get_header_attr(self, attr):
        """Return an attribute of the simulation header
        It has alredy been loaded in the constructor"""
        return self.header[attr]
    
    def _load_block_header(self, part_type, blockname):
        """Get the header of a block
        e.g. PART_272/0/Position/attr-v2
        """
        block_path = self._get_block_path(part_type, blockname)
        header = {'dtype': None,
                    'nmemb': None,
                    'nblobs': None,
                    'blob_names':[],
                    'blob_part_start':None,
                    'blob_part_end':None}
        # We only used the 0 , 2nd column of the header
        # But, get others here just in case
        # Column 0: blob name in hexadecimal
        # Column 1: number of particles in the blob
        # Column 2: the checksum of the blob
        # Column 3: A mysterious number Yu uses here:
        # https://github.com/MP-Gadget/bigfile/blob/master/src/bigfile.c#L586
        data = [] # Store all 3 columns as integers
        header_path = os.path.join(block_path, 'header')
        assert os.path.exists(header_path)
        with open(header_path) as f:
            lines = f.readlines()
        # Get the dtype, dimension, and numbed of blobs
        header['dtype'] = np.dtype(lines[0].split(':')[1].strip())
        header['nmemb'] = int(lines[1].split(':')[1].strip())
        header['nblobs'] = int(lines[2].split(':')[1].strip())
        # Initialize an empty list to store the data
        data = []
        # Skip the first three header lines
        data_lines = lines[3:]
        # Loop through each data line
        for line in data_lines:
            # Remove leading/trailing whitespace
            line = line.strip()
            # Split the line at the first colon to separate the index
            idx_hex, rest = line.split(':', 1)
            header['blob_names'].append(idx_hex)
            # Convert the hexadecimal index to an integer
            idx_int = int(idx_hex.strip(), 16)
            # Split the rest of the line at ' : ' to get the data columns
            values = rest.strip().split(' : ')
            # Convert the string values to integers
            values = [int(v.strip()) for v in values]
            # Combine the index and the values into one list
            row = [idx_int] + values
            # Append the row to the data list
            data.append(row)
        # Convert the data list to a NumPy array
        data = np.array(data)
        assert data.shape == (header['nblobs'], 4), f"header shape is {data.shape}"
        # Find the start and end particle numbers for each blob
        header['blob_part_end'] = np.cumsum(data[:,1])
        blob_part_start = np.concatenate(([0], header['blob_part_end'][:-1]))
        header['blob_part_start'] = blob_part_start

        return header


    def _get_block_header_attr(self, part_type, blockname): 
        """Get the header of a block which is already loaded"""
        if os.path.join(str(part_type), blockname) in self.block_headers.keys():
            return self.block_headers[os.path.join(str(part_type), blockname)]
        else:
            # load the block header on root rank and broadcast it
            header = None
            if self.rank == 0:
                header = self._load_block_header(part_type, blockname)
            if self.comm is not None:
                header = self.comm.bcast(header, root=0)
            self.block_headers[os.path.join(str(part_type), blockname)] = header
            return header

    def _get_block_path(self, part_type, blockname):
        """Get the path to a block"""
        if blockname in self.hdf_to_bigfile_map.keys():
            blockname = self.hdf_to_bigfile_map[blockname]
        block_path = os.path.join(self.part_dir, str(part_type), blockname)
        return block_path

    def _get_blobs_for_rank(self, segment, part_type, blockname):
        """Find the proper blobs having the particles this segment of the rank needs
        Parameters:
            segment (int): The segment number
            part_type (int): The particle type
            blockname (str): The block name, e.g. 'Position'
        Returns:
            blobs (np.ndarray): The blob indices
            blob_paths (list): The paths to the blobs
            start_blob (np.ndarray): The start particle number for each blob, indexing from 0
            for each blob
            end_blob (np.ndarray): The end particle number for each blob, indexing from 0
            for each blob
            """
        # Get the header of the block
        header = self._get_block_header_attr(part_type, blockname)
        (start, end) = self._segment_to_partlist(part_type = part_type, segment=segment)
        # Get the list of files for this rank
        first = np.where(header['blob_part_start'] <= start)[0][-1]
        last = np.where(header['blob_part_end'] >= end)[0][0]
        # The blob ids containing the particles
        blobs = np.arange(first, last+1)
        
        start_blob = start - header['blob_part_start'][blobs]
        start_blob[start_blob < 0] = 0  # Start from the beginning of the blob
        end_blob = end - header['blob_part_start'][blobs]

        ind = np.where(header['blob_part_end'][blobs] < end)[0]
        end_blob[ind] = (header['blob_part_end'] - header['blob_part_start'])[blobs][ind]
        
        blob_paths = [os.path.join(self._get_block_path(part_type, blockname), header['blob_names'][b]) for b in blobs]
        return blobs, blob_paths, start_blob, end_blob

    def get_data(self, part_type, blockname, segment):
        """Get the data for a particular block, specified
        using either the BigFile names or the HDF5 names.
        Segment: which data segment to load."""
        if blockname in self.hdf_to_bigfile_map.keys():
            blockname = self.hdf_to_bigfile_map[blockname]
        
        blob_ids, blob_paths, start_blobs, end_blobs = self._get_blobs_for_rank(segment, part_type, blockname)
        dtype = self._get_block_header_attr(part_type, blockname)['dtype']
        # The dimension of the block, i.e. (# Parts, nmembs)
        nmembs = self._get_block_header_attr(part_type, blockname)['nmemb']
        # particles this rank is responsible for
        data = np.array([])
        for i, bl_id in enumerate(blob_ids):
            # Offset in the blob file to start reading from
            offset = start_blobs[i] * nmembs * dtype.itemsize
            # The buffer to read the particles from this blob
            buffer = np.empty( (end_blobs[i] - start_blobs[i])*nmembs, dtype=dtype)
            # MPI-IO read the particles from the file
            # Note: Since each file won't be openned by many ranks, 
            # i.e. `comm.size < num blob files`, no need for collective 
            # synchronization among COMM_WORLD and each rank opens the 
            # blobs with COMM_SELF
            f_handle = self.MPI.File.Open(self.MPI.COMM_SELF, blob_paths[i], self.MPI.MODE_RDONLY)
            f_handle.Read_at(offset, buffer)
            data = np.append(data, np.frombuffer(buffer, dtype= dtype))
            f_handle.Close()
            del f_handle
            del buffer
        
        data = data.reshape((data.size//nmembs, nmembs))
        return data
    
    def get_n_segments(self, part_type, chunk_size = 512.**3):
        """Distribute particles among ranks. Also, break the load on each rank into segments containing
        chunk_size particles and return number of these segments for each rank""" 
        #Store number of particles on each rank in an array
        self.parts_rank = ((self.get_npart()[part_type]//self.size)*np.ones(shape=(self.size,))).astype(int)
        remainder = int(self.get_npart()[part_type]%self.size)
        # Some ranks would get one more particle
        if remainder !=0 :
            self.parts_rank[0:remainder] += 1
        # Have an extra chunk if the number of particles is not a multiple of chunk_size
        nseg = int(np.ceil(np.max([1,self.parts_rank[self.rank]/chunk_size])))
        return nseg, chunk_size
    
    def get_blocklen(self, part_type, blockname, segment):
        """Get the length of a block"""
        if blockname in self.hdf_to_bigfile_map.keys():
            blockname = self.hdf_to_bigfile_map[blockname]
        try:
            (start, end) = self._segment_to_partlist(part_type=part_type, segment=segment)
            if end is not None:
                return end - start
            #Last segment has no end
            return self._f_handle[str(part_type)+"/"+blockname].size - start

        except bigfile.BigFileError:
            raise KeyError("Not found:"+str(part_type)+"/"+blockname)

    def _segment_to_partlist(self, part_type, segment):
        """Get the first and last particle in a segment.
        Parameters:
            part_type (int): The particle type
            segment (int): The segment number
        Returns:
            (int, int): The first and last particle in the
                        segment for this rank.
        """
        if segment < 0:
            return (0, None)
        # Get number of segments on this rank which we break particles into
        n_segments, chunk_size = self.get_n_segments(part_type)
        if self.rank ==0 :
            first_part_rank = 0
        else:
            # First particle for this rank depends on particles taken by previous ranks
            first_part_rank = np.sum(self.parts_rank[0:self.rank])
        start_seg = first_part_rank+int(chunk_size)*segment
        end_seg = start_seg + int(chunk_size)
        # The last chunk may have fewer particles than chunk_size
        if end_seg-first_part_rank > self.parts_rank[self.rank]:
            end_seg = self.parts_rank[self.rank]+first_part_rank
        return (start_seg, end_seg)

    def get_kernel(self):
        """Get the integer corresponding to a density kernel for each particle.
           SPH cubic and quintic splines are implemented. Values are:
               0 - Top hat kernel (for Arepo)
               1 - SPH cubic spline kernel (for Gadget).
               2 - A partial reconstruction of the Voronoi mesh, for Arepo.
               3 - SPH quintic spline kernel (for modern SPH Gadget)
            The types defined in allvars.h of MP-Gadget are:
                DENSITY_KERNEL_CUBIC_SPLINE = 1,
                DENSITY_KERNEL_QUINTIC_SPLINE = 2,
                DENSITY_KERNEL_QUARTIC_SPLINE = 4,
        """
        try:
            kernel = self.get_header_attr("DensityKernel")
            if kernel == 2:
                return 3
        except KeyError:
            kernel = 1
        #Other types are not yet supported.
        assert kernel == 1 or kernel == 3
        return kernel

    def get_peculiar_velocity(self, part_type, segment):
        """Get the peculiar velocity in internal units. Converts out the various Gadget comoving a factors."""
        vel = self.get_data(part_type, "Velocity", segment = segment)
        pecvel = self.get_header_attr("UsePeculiarVelocity")
        if not pecvel:
            atime = self.get_header_attr("Time")
            vel /= atime
        return vel
