# -*- coding: utf-8 -*-
"""A small module for computing the smoothing length of a Gadget/Arepo simulation."""

import numpy as np

def get_smooth_length(bar):
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
        radius = np.cbrt(bar["Volume"], dtype=np.float32)
    except KeyError:
        #If we don't have a Volume array we are gadget, and
        #the SmoothingLength array is actually the smoothing length.
        #There is a different kernel definition, as in gadget the kernel goes from 0 to 2,
        #whereas I put it between zero and 1.
        radius=np.array(bar["SmoothingLength"],dtype=np.float32)/2
    except AttributeError:
        #This is for really old numpys without cbrts
        radius = np.power(bar["Volume"], 1./3, dtype=np.float32)

    return radius
