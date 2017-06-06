# -*- coding: utf-8 -*-
"""Module to compute flux statistics from spectra:
the power spectrum, the pdf and to normalise to a mean tau.
Useful for lyman alpha forest work."""

import math
import numpy as np
from scipy.optimize import brentq

from ._spectra_priv import _rescale_mean_flux

def obs_mean_tau(redshift):
    """The mean flux from 0711.1862: is (0.0023±0.0007) (1+z)^(3.65±0.21)
    Todo: check for updated values."""
    return 0.0023*(1.0+redshift)**3.65

def mean_flux(tau, mean_flux_desired, tol = 1e-5):
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
    return _rescale_mean_flux(tau, mean_flux_desired, np.size(tau), tol)

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
    #Normalise the FFT: equivalent to using norm = ortho for recent numpy versions
    power /= np.shape(inarray)[axis]
    return power

def flux_power(tau, vmax, mean_flux_desired=None):
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
    else:
        mean_flux_desired = np.mean(np.exp(-tau))
    (_, npix) = np.shape(tau)
    dflux=np.exp(-scale*tau)/mean_flux_desired - 1.
    # Calculate flux power for each spectrum in turn
    flux_power_perspectra = _powerspectrum(dflux, axis=1)
    #Take the mean and convert units.
    mean_flux_power = vmax*np.mean(flux_power_perspectra, axis=0)
    assert np.shape(mean_flux_power) == (npix//2+1,)
    kf = _flux_power_bins(vmax, npix)
    return kf,mean_flux_power

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
