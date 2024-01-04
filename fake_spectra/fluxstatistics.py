# -*- coding: utf-8 -*-
"""Module to compute flux statistics from spectra:
the power spectrum, the pdf and to normalise to a mean tau.
Useful for lyman alpha forest work."""

import math
import numpy as np
from datetime import datetime

from nbodykit.lab import FFTPower
from nbodykit.source.catalog import ArrayCatalog
from nbodykit import setup_logging
from ._spectra_priv import _rescale_mean_flux

def obs_mean_tau(redshift):
    """The mean flux from 0711.1862: effective optical depth is (0.0023±0.0007) (1+z)^(3.65±0.21)
    Todo: check for updated values."""
    return 0.0023*(1.0+redshift)**3.65

def mean_flux(tau, mean_flux_desired, tol = 1e-5, thresh=1e30):
    """Scale the optical depths by a constant value until we get the observed mean flux.
    ie, we want F_obs = bar{F} = < e^-tau >
    Solves iteratively using Newton-Raphson.
    This is safe because the exponential function is so well-behaved.
    Arguments:
        tau - optical depths to scale
        mean_flux_desired - mean flux desired
        tol - tolerance within which to hit mean flux
    returns:
        scaling factor for tau"""
    if np.size(tau) == 0:
        return 0
    return _rescale_mean_flux(tau.astype(np.float64), mean_flux_desired, np.size(tau), tol, thresh)

def flux_pdf(tau, nbins=20, mean_flux_desired=None):
    """Compute the flux pdf, a normalised histogram of the flux, exp(-tau)"""
    scale = 1.
    if mean_flux_desired is not None:
        scale = mean_flux(tau, mean_flux_desired)
    flux = np.exp(-scale * tau)
    bins = np.arange(nbins+1)/(1.*nbins)
    (fpdf, _) = np.histogram(flux, bins=bins,density=True)
    cbins = (bins[1:] + bins[:-1])/2.
    return cbins, fpdf

def _powerspectrum(inarray, axis=-1):
    """Compute the power spectrum of the input using np.fft"""
    rfftd = np.fft.rfft(inarray, axis=axis)
    # Want P(k)= F(k).re*F(k).re+F(k).im*F(k).im
    power = np.abs(rfftd)**2
    #Normalise the FFT so it is independent of input size.
    power /= np.shape(inarray)[axis]**2
    return power

def _window_function(k, *, R, dv):
    """The window function corresponding to the spectra response of the spectrograph.
    R is the spectrograph resolution.
    dv is the pixel width of the spectrograph.
    Default values for BOSS are:
        dv = 69, R = 60 at 5000 A and R = 80 at 4300 A."""
    #FWHM of a Gaussian is 2 \sqrt(2 ln 2) sigma
    sigma = R/(2*np.sqrt(2*np.log(2)))
    return np.exp(-0.5 * (k * sigma)**2) * np.sinc(k * dv/2/math.pi)

def flux_power(tau, vmax, spec_res = 8, mean_flux_desired=None, window=False):
    """Get the power spectrum of (variations in) the flux along the line of sight.
        This is: P_F(k_F) = <d_F d_F>
                 d_F = e^-tau / mean(e^-tau) - 1
        If mean_flux_desired is set, the spectral optical depths will be rescaled
        to match the desired mean flux.
        We compute the power spectrum along each sightline and then average the result.
        Arguments:
            tau - optical depths. Shape is (NumLos, npix)
            mean_flux_desired - Mean flux to rescale to.
	    vmax - velocity scale corresponding to maximal length of the sightline.
        Returns:
            flux_power - flux power spectrum in km/s. Shape is (npix)
            bins - the frequency space bins of the power spectrum, in s/km.
    """
    scale = 1.
    if mean_flux_desired is not None:
        scale = mean_flux(tau, mean_flux_desired)
        #print("rescaled: ",scale,"frac: ",np.sum(tau>1)/np.sum(tau>0))
    else:
        mean_flux_desired = np.mean(np.exp(-tau))
    (nspec, npix) = np.shape(tau)
    mean_flux_power = np.zeros(npix//2+1, dtype=tau.dtype)
    # compute in batches, purely for computational efficiency
    for i in range(10):
        end = min((i+1)*nspec//10, nspec)
        dflux=np.exp(-scale*tau[i*nspec//10:end])/mean_flux_desired - 1.
        # Calculate flux power for each spectrum in turn
        flux_power_perspectra = _powerspectrum(dflux, axis=1)
        #Take the mean and convert units.
        mean_flux_power += vmax*np.sum(flux_power_perspectra, axis=0)
    mean_flux_power/= nspec
    assert np.shape(mean_flux_power) == (npix//2+1,)
    kf = _flux_power_bins(vmax, npix)
    #Divide out the window function
    if window and spec_res > 0:
        mean_flux_power /= _window_function(kf, R=spec_res, dv=vmax/npix)**2
    return kf,mean_flux_power

def _3d_powerspectrum(dflux_mesh, boxsize, los, dk=None, Nmu=6):
    """Compute the 3D power spectrum of the input using nbodykit
    Parameters:
    dfux - 3D array of flux variations
    boxsize - size of the box in units of interest (eg, comoving cMpc/h), 
                the units of the 3d power spectrum, i.e. P(k,mu), will be in these units
    los - line of sight direction, i.e. [0,0,1] for z-axis
    dk - bin width in k
    Nmu - number of mu bins
    Returns:
    power - a dictionary with the p(k,mu) and the k and mu bins, keys:['power','k','mu']
    """
    power = FFTPower(dflux_mesh, BoxSize=boxsize, mode='2d', los= los, dk=dk, Nmu=Nmu,
                     save_3d_power=False)
    return power.power

def flux_power_3d(tau, boxsize, mean_flux_desired=None, dk=None, Nmu=6):
    """Get the power spectrum of (variations in) the flux in 3D which is binned in (k,mu).
        This is: P_3D(k) = <d_F d_F>
                 d_F = e^-tau / mean(e^-tau) - 1
                 Then we bin P_3D(k) in k and mu.
        If mean_flux_desired is set, the spectral optical depths will be rescaled
        to match the desired mean flux.
        We compute the power spectrum along each sightline and then average the result.
        Arguments:
            tau - optical depths. Shape is (NumLos, npix)
            mean_flux_desired - Mean flux to rescale to.
        boxsize - size of the box in units of interest (eg, comoving cMpc/h), 
                the units of the 3d power spectrum, i.e. P(k,mu), will be in these units
        Returns:
            k, mu - the k and mu bins of the power spectrum
            flux_power - flux power spectrum in `boxsize` units
            Note: The first row corresponds to k=0, so you can remove it later
    """
    scale = 1.
    if mean_flux_desired is not None:
        scale = mean_flux(tau, mean_flux_desired)
        #print("rescaled: ",scale,"frac: ",np.sum(tau>1)/np.sum(tau>0))
    else:
        mean_flux_desired = np.mean(np.exp(-tau))
    (nspec, npix) = np.shape(tau)
    nt = np.sqrt(nspec/3).astype(int)
    x, y, z = np.meshgrid(np.arange(nt), np.arange(nt), np.arange(npix), indexing='ij')
    x = x*boxsize/nt
    y = y*boxsize/nt
    z = z*boxsize/npix
    coords = np.vstack((x.ravel(), y.ravel(), z.ravel())).T
    for i in range(3):
        print(f'Interpolating the spectra along the perp direction | {datetime.now()}', flush=True)
        end = min((i+1)*nspec//3, nspec)
        #Interpolate onto a regular grid
        setup_logging()
        cat = ArrayCatalog({'Position': coords, 'Weight': tau[i*nspec//3:end].ravel()})
        mesh = cat.to_mesh(Nmesh=[nt, nt, nt], BoxSize=boxsize, weight='Weight', resampler='cic')
        get_df = lambda x, v: np.exp(-scale*v)/mean_flux_desired - 1
        mesh = mesh.apply(get_df, mode='real', kind='index')
        # Calculate flux power for the interpolated flux field
        print(f'Calculating the 3D power spectrum for axis {i} | {datetime.now()}', flush=True)
        los = [0, 0, 0]
        los[i] = 1
        power = _3d_powerspectrum(mesh, boxsize=boxsize, los=los, dk=dk, Nmu=Nmu)
        #The units of the P(k,mu) is same as the `boxsize` argument
        if i==0:
            mean_flux_power = power['power']
        else:
            mean_flux_power += power['power']
    # Avergaing the 3D power spectrum obtained with the spectra along the 3 axes
    mean_flux_power/= 3
    # nobodykit calculates the power in boxsize unit
    k = power['k']
    mu = power['mu']
    # We do not do any window correction along los or the transverse directions
    # We have seen this effect been marginal for the 1D power spectrum becase
    # the spatial resolution is much larger than the observations
    return k, mu, mean_flux_power

def _flux_power_bins(vmax, npix):
    """
        Generate k bins for the flux power spectrum by converting the natural
        (ie, fractions of the total spectrum) units output by the flux power spectrum
        routine into physical km/s, accounting for Fourier convention.
        Arguments:
            vmax - the length of a spectrum in km/s and the conversion factor from comoving kpc is:
                H(z) * a / h / 1000
                defined in spectra.py:115
            nbins - number of bins of *input spectrum* - not the fourier output!
        Returns: bin center in s/km
    """
    #Get the frequency component
    kf = np.fft.rfftfreq(npix)
    #Units:
    #The largest frequency scale is the velocity scale of the box,
    #not 1/nbins as rfftfreq gives.
    #Adjust Fourier convention.
    kf *= 2.0*math.pi * npix/vmax
    return kf
