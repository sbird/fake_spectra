# -*- coding: utf-8 -*-
"""Split off module to load halo catalogues and export a list of mass and positions"""

import numpy as np
from . import subfindhdf
try:
    xrange(1)
except NameError:
    xrange = range

#Internal gadget mass unit: 1e10 M_sun/h in g/h
UnitMass_in_g=1.989e43
#1 M_sun in g
SolarMass_in_g=1.989e33


def is_masked(halo,sub_mass,sub_cofm, sub_radii):
    """Find out whether a halo is a mere satellite and if so mask it"""
    near=np.where(np.all((np.abs(sub_cofm[:,:]-sub_cofm[halo,:]) < sub_radii[halo]),axis=1))
    #If there is a larger halo nearby, mask this halo
    return np.size(np.where(sub_mass[near] > sub_mass[halo])) == 0

def find_all_halos(num, base, min_mass):
    """Get a halo catalogue and return its members, filtering out those with masses below min_mass.
    Select halos via their M_200 mass, defined in terms of the critical density.
    Arguments:
        num - snapnumber
        base - simulation directory
        min_mass - minimum mass of halos to use
    Returns:
        ind - list of halo indices used
        sub_mass - halo masses in M_sun /h
        sub_cofm - halo positions
        sub_radii - R_Crit200 for halo radii"""
    subs=subfindhdf.SubFindHDF5(base, num)
    #Get list of halos resolved, using a mass cut; cuts off at about 2e9 for 512**3 particles.
    ind=np.where(subs.get_grp("Group_M_Crit200") > min_mass)
    #Store the indices of the halos we are using
    #Get particle center of mass, use group catalogue.
    sub_cofm=np.array(subs.get_grp("GroupPos")[ind])
    #halo masses in M_sun/h: use M_200
    sub_mass=np.array(subs.get_grp("Group_M_Crit200")[ind])*UnitMass_in_g/SolarMass_in_g
    #r200 in kpc/h (comoving).
    sub_radii = np.array(subs.get_grp("Group_R_Crit200")[ind])
    del subs

    return (ind, sub_mass,sub_cofm,sub_radii)

def find_wanted_halos(num, base, min_mass, dist=1):
    """When handed a halo catalogue, remove from it the halos that are within dist virial radii of other, larger halos.
    Select halos via their M_200 mass, defined in terms of the critical density.
    Arguments:
        num - snapnumber
        base - simulation directory
        min_mass - minimum mass of halos to use
        dist - Factor to multiply the virial radius by
    Returns:
        ind - list of halo indices used
        sub_mass - halo masses in 1e10 M_sun /h
        sub_cofm - halo positions
        sub_radii - dist*R_Crit200 for halo radii"""

    (ind, sub_mass,sub_cofm,sub_radii) = find_all_halos(num, base, min_mass)
    sub_radii*=dist
    #For each halo
    ind2=np.where([is_masked(ii,sub_mass,sub_cofm,sub_radii) for ii in xrange(0,np.size(sub_mass))])
    ind=(np.ravel(ind)[ind2],)
    sub_mass=sub_mass[ind2]
    sub_cofm=sub_cofm[ind2]
    sub_radii=sub_radii[ind2]
    return (ind, sub_mass,sub_cofm,sub_radii)

