# -*- coding: utf-8 -*-
"""A small module for computing the smoothing length of a particle simulation.
(Non-trivial in Arepo)"""

from math import pi
import numpy as np
import numexpr as ne

#To get these factors as floats for numexpr
scal = np.float32(3./4./pi)
poww = np.float32(1./3.)


def get_smooth_length(bar):
    """Figures out if the particles are from AREPO or GADGET
    and computes the smoothing length.
    Note the Volume array in HDF5 is comoving and this returns a comoving smoothing length
    If we are Arepo, this smoothing length is  cell radius, where
    cell volume = 4/3 π (cell radius) **3 and cell volume = mass / density
    Arguments:
        Baryon particles from a simulation
    Returns:
        Array of smoothing lengths in code units.
    """
    #Are we arepo? If we are a modern version we should have this array.
    if 'Volume' in bar.keys():
        volume=np.array(bar["Volume"],dtype=np.float32)
        radius = ne.evaluate("(scal*volume)**poww")
    elif 'Number of faces of cell' in bar.keys():
        rho=np.array(bar["Density"])
        mass=np.array(bar["Masses"])
        radius = ne.evaluate("(scal*mass/rho)**poww")
    else:
        #If we are gadget, the SmoothingLength array is actually the smoothing length.
        radius=np.array(bar["SmoothingLength"])
    return radius
