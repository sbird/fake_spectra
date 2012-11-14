import numpy as np
import hsml

from _spectra_priv import _SPH_Interpolate

def SPH_Interpolate(data, los_table, nbins, box):
    """Interpolate particles to lines of sight, calculating density, temperature and velocity
    of various species (TODO: only hydrogen now) along the line of sight.

    This is a wrapper which calls the C function.
    Arguments:
    	data - HDF5 dataset from snapshot. Use f["PartType0"]
	los_table - table of los positions. should have member arrays x, y, z and axis.
	nbins - number of bins in each spectrum
	box - box size
    """
    pos = np.array(data["Coordinates"],dtype=np.float32)
    vel = np.array(data["Velocities"],dtype=np.float32)
    mass = np.array(data["Masses"],dtype=np.float32)
    u = np.array(data["InternalEnergy"],dtype=np.float32)
    nh0 = np.array(data["NeutralHydrogenAbundance"],dtype=np.float32)
    ne = np.array(data["ElectronAbundance"],dtype=np.float32)
    hh = np.array(hsml.get_smooth_length(data),dtype=np.float32)
    xx=np.array(los_table.xx, dtype=np.float32)
    yy=np.array(los_table.yy, dtype=np.float32)
    zz=np.array(los_table.zz, dtype=np.float32)
    axis=np.array(los_table.axis, dtype=np.int32)
    return _SPH_Interpolate(nbins, box, pos, vel, mass, u, nh0, ne, hh, axis, xx, yy, zz)
