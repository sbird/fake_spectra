# -*- coding: utf-8 -*-
"""
Python module for generating fake total emission in a magnitude band along a sightline.

This uses the Arepo/Illustris output GFM_Photometrics to get photometric band data,
which may or may not be accurate.
"""

from __future__ import print_function
import math
import os.path as path
import shutil
import h5py
import numpy as np

from . import spectra as ss

def maginJy(mag, band):
    """Convert a magnitude to flux in Jansky, according to wikipedia's table"""
    bandfluxes = {'U':1810, 'B':4260, 'V':3640,'K':670,'g':3730,'r':4490,'i':4760 ,'z':4810}
    return 10**(mag/(-2.5))*bandfluxes[band]

def apparentflux(DL):
    """Convert flux from absolute magnitudes (flux at 10 pc distance) to apparent flux in Jy.
    DL is luminosity distance in Mpc"""
    return (10/(DL*1e6))**2

def distance(arcsec, redshift, hubble, OmegaM):
    """Find the size of something in comoving kpc/h from the size on the sky in arcseconds.
    """
    #First arcsec to radians
    #2 pi radians -> degrees -> arcminute -> arcsecond
    rad = 2*math.pi/360./60./60. * arcsec
    #Then to physical kpc
    atime = 1./(1+redshift)
    (_, DA, _) = calculator(hubble*100, OmegaM, redshift)
    size = DA * rad * 1000
    #Comoving kpc/h
    size = size /( atime/ hubble)
    return size

class EmissionSpectra(ss.Spectra):
    """Class to compute the emission from stars in B band around the DLA spectrum"""
    stellar = {}

    def _read_stellar_data(self,fn, band, hhmult=10.):
        """Read the particle data for a single interpolation"""
        bands = {'U':0, 'B':1, 'V':2,'K':3,'g':4,'r':5,'i':6,'z':7}
        nband = bands[band]
        pos = self.snapshot_set.get_data(4,"Position", segment = fn).astype(np.float32)
        #Set each stellar radius to the pixel size
        hh = hhmult*np.ones(np.shape(pos)[0], dtype=np.float32)
        #Find particles we care about
        ind = self.particles_near_lines(pos, hh,self.axis,self.cofm)
        #print np.size(ind)
        #Do nothing if there aren't any, and return a suitably shaped zero array
        if np.size(ind) == 0:
            raise ValueError("No stars")
        pos = pos[ind,:]
        hh = hh[ind]
        #Find the magnitude of stars in this band
        emflux = maginJy(self.snapshot_set.get_data(4,"GFM_StellarPhotometrics", segment = fn).astype(np.float32)[ind][:,nband],band)
        fluxx = np.array([ np.sum(emflux[self.particles_near_lines(pos, hh,np.array([ax,]),np.array([cofm,]))]) for (ax, cofm) in zip(self.axis, self.cofm)])
        #print np.sum(emflux)
        return fluxx
        #return (pos, emflux, hh)

    def get_emflux(self, band, pixelsz=1):
        """
        Get the density weighted flux in each pixel for a given species.
        band: rest-frame optical band observed in
        pixelsz: Angular size of the pixels in arcseconds
        """
        #Mapping from bandname to number
        dist = distance(pixelsz, 1./self.atime-1, self.hubble, self.OmegaM)
        try:
            self._really_load_array((band,dist), self.stellar, "stellar")
            emflux = self.stellar[(band,dist)]
        except KeyError:
            emflux = np.zeros(self.NumLos,dtype=np.float32)
            for fn in self.snapshot_set.get_n_segments():
                try:
                    emflux +=  self._read_stellar_data(fn, band,dist)
                except ValueError:
                    pass
            self.stellar[(band,dist)] = emflux
        (_,_,DL) = calculator(self.hubble*100, self.OmegaM, 1./self.atime-1)
        emflux *= apparentflux(DL)
        return emflux

    def save_file(self):
        """
        Saves spectra to a file, because they are slow to generate.
        File is by default to be $snap_dir/snapdir_$snapnum/spectra.hdf5.
        """
        #We should make sure we have loaded all lazy-loaded things first.
        self._load_all_multihash(self.stellar, "stellar")
        self._load_all_multihash(self.tau_obs, "tau_obs")
        self._load_all_multihash(self.tau, "tau")
        self._load_all_multihash(self.colden, "colden")
        try:
            self._load_all_multihash(self.colden, "velocity")
        except IOError:
            pass
        try:
            if path.exists(self.savefile):
                shutil.move(self.savefile,self.savefile+".backup")
            f=h5py.File(self.savefile,'w')
        except IOError:
            try:
                f=h5py.File(self.savefile,'w')
            except IOError:
                raise IOError("Could not open ",self.savefile," for writing")
        grp_grid = f.create_group("stellar")
        self._save_multihash(self.stellar, grp_grid)
        self._save_file(f)


def calculator(H0, Omega_M, zz):
    """Compute luminosity distance for a given cosmology. Assumes flatness.
        Freely adapted from James Schombert's python version of Ned Wright's cosmology calculator.
    Inputs:
        H0 - Hubble constant in km/s/Mpc
        Omega_M - Omega_matter
        zz - redshift to compute distances to
    Returns:
        (Comoving distance, angular distancem luminosity distance) (all in physical Mpc)"""

    light = 299792.458 # speed of light in km/sec
    h = H0/100.
    WR = 4.165E-5/(h*h)   # includes 3 massless neutrino species, T0 = 2.72528
    #Assume flat
    WV = 1- Omega_M-WR
    #scale factor to compute distance to
    az = 1.0/(1.+zz)
    n=1000         # number of points in integrals

    # do integral over a=1/(1+z) from az to 1 in n steps, midpoint rule
    a = np.logspace(np.log10(az), 0, n)
    a2H = a*a*np.sqrt(Omega_M/a**3+WR/(a**4)+WV)
    #Comoving distance
    DCMR = np.trapz(1./a2H, a)
    #In Mpc
    DC_Mpc = (light/H0) * DCMR
    # angular size distance In Mpc
    DA_Mpc = (light/H0)*az*DCMR
    #Luminosity distance in Mpc
    DL_Mpc = DA_Mpc/(az*az)

    #print 'The comoving radial distance is %1.1f' % DC_Mpc + ' Mpc'
    #print 'The angular size distance D_A is ' + '%1.1f' % DA_Mpc + ' Mpc'
    #print 'The luminosity distance D_L is ' + '%1.1f' % DL_Mpc + ' Mpc'
    return (DC_Mpc, DA_Mpc, DL_Mpc)
