# -*- coding: utf-8 -*-
"""Module to compute flux statistics from spectra:
the power spectrum, the pdf and to normalise to a mean tau.
Mostly useful for lyman alpha studies."""

import numpy as np
import math
from scipy.optimize import brentq

def flux_power(tau, redshift, rescale=True, statistic = 'power',vmax=None):
    """Compute a flux statistic, potentially rescaling to the observed mean flux.
        Arguments:
            tau - optical depths
            redshift - redshift at which to use the observed mean flux.
            statistic - If power computes the flux power spectrum. If pdf computes the flux pdf.
            rescale - Shall I rescale the optical depths?
        Returns:
            Computed flux power or pdf.
    """
    obs_flux = obs_mean_flux(redshift)
    if rescale:
        scale = mean_flux(tau, obs_flux)
    else:
        scale = 1.
    if statistic == 'power':
        if vmax == None:
            raise ValueError("You must specify vmax to get the power spectrum")
        pk = _calc_power(tau, obs_flux, scale)*vmax
        bins = _flux_power_bins(vmax, np.shape(tau)[1])
        return (pk, bins)
    elif statistic == 'pdf':
        stat = _calc_pdf(tau, scale=scale)
        return stat
    else:
        raise NotImplementedError("Only flux power and pdf implemented")


def obs_mean_flux(redshift):
    """The mean flux from 0711.1862: is (0.0023±0.0007) (1+z)^(3.65±0.21)
    Todo: check for updated values."""
    return 0.0023*(1.0+redshift)**3.65

def _mean_flux_kernel(scale, tau, obs_flux):
    """Kernel for scipy optimisation routine below"""
    return np.mean(np.exp(-scale*tau)) - obs_flux

def mean_flux(tau, obs_flux, tol = 0.05):
    """Scale the optical depths by a constant value until we get the observed mean flux.
    Solves implicitly using Newton-Raphson.
    This is safe because the exponential function is so well-behaved.
    Arguments:
        tau - optical depths to scale
        obs_flux - mean flux desired
        tol - tolerance within which to hit mean flux
    returns:
        scaling factor for tau"""
    #Find the initial interval by simple doubling
    scale = 10
    while scale < 100:
        if _mean_flux_kernel(scale, tau, obs_flux) < 0:
            break
        scale*=2
    newscale = brentq(_mean_flux_kernel, 0, scale,args=(tau, obs_flux), xtol=tol)
    return newscale

def _calc_pdf(tau, pbins=20, scale=1.):
    """Compute the flux pdf, a normalised histogram of the flux, exp(-tau)"""
    return np.histogram(np.exp(-scale*tau), bins = pbins, range=(0,1), density=True)

def _powerspectrum(inarray, axis=-1):
    """Compute the power spectrum of the input using np.fft"""
    rfftd = np.fft.rfft(inarray, axis=axis)
    # Want P(k)= F(k).re*F(k).re+F(k).im*F(k).im
    power = np.abs(rfftd)**2
    return power

def _calc_power(tau, tau_eff, scale=1.):
    """
        Compute the flux power spectrum along the line of sight.
        This is really the power spectrum of delta_flux = exp(-tau) / exp(-tau_eff) - 1
        We compute the power spectrum along the line of sight to each quasar and then average the result.
        Arguments:
            tau - optical depths. Shape is (NumLos, npix)
            tau_eff - Mean flux.
            nbins - Number of bins for the output power spectrum.
        Returns:
            flux_power - flux power spectrum. Shape is (npix)
            bins - the frequency space bins of the power spectrum, in s/km.
    """
    (_, npix) = np.shape(tau)
    dflux=np.exp(-scale*tau)/np.exp(-tau_eff)-1.
    # Calculate flux power for each spectrum in turn
    flux_power_perspectra = _powerspectrum(dflux, axis=1)
    #Take the mean
    mean_flux_power = np.mean(flux_power_perspectra, axis=0)
    assert np.shape(mean_flux_power) == (npix//2+1,)
    return mean_flux_power

def _flux_power_bins(vmax, nbins):
    """
        Generate k bins for th eflux power spectrum by converting the natural
        (ie, fractions of the total spectrum) units output by the flux power spectrum
        routine into physical km/s, accounting for Fourier convention.
        Arguments:
            vmax - the length of a spectrum in km/s and the conversion factor from comoving kpc is:
                H(z) * a / h / 1000
                defined in spectra.py:115
            nbins - number of bins of *input spectrum* - not the fourier output!
        Returns: bin center in s/km
    """
    kk = (np.arange(0,nbins//2+1)+0.5)/vmax*2*math.pi
    return kk
