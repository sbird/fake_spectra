# -*- coding: utf-8 -*-
"""Module to compute flux statistics from spectra: 
the power spectrum, the pdf and to normalise to a mean tau.
Mostly useful for lyman alpha studies."""

import spectra as ss
import numpy as np

#class Statistics(ss.Spectra):

def obs_mean_flux(redshift):
    """The mean flux from 0711.1862: is (0.0023±0.0007) (1+z)^(3.65±0.21)
    Todo: check for updated values."""
    return 0.0023*(1.0+redshift)**3.65
    
def mean_flux(tau, obs_flux, tol = 0.05):
    """Scale the optical depths by a constant value until we get the observed mean flux.
    Solves implicitly.
    Arguments:
        tau - optical depths to scale
        obs_flux - mean flux desired
        tol - tolerance within which to hit mean flux
    returns: 
        scaling factor for tau"""
    newscale=1.
    scale = 0.
    while abs(newscale-scale) > abs(tol*newscale)+1e-6:
        scale=newscale
        flux = np.exp(-scale*tau)
        mean_flux = np.mean(flux)
        tau_mean_flux = np.mean(tau*flux)
        newscale=scale+(mean_flux-obs_flux)/tau_mean_flux;
    return newscale

def calc_pdf(tau, pbins=20):
    """Compute the flux pdf, a normalised histogram of the flux, exp(-tau)"""
    return np.histogram(np.exp(-tau), bins = pbins, range=(0,1), density=True)

def powerspectrum(inarray, axis=-1):
    """Compute the power spectrum of the input using np.fft"""
    rfftd = np.fft.rfft(inarray, axis=axis)
    # Want P(k)= F(k).re*F(k).re+F(k).im*F(k).im
    power = np.abs(rfftd)**2
    return power

def calc_power(tau, tau_eff):
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
    """
    (_, npix) = np.shape(tau)
    dflux=np.exp(-tau)/np.exp(-tau_eff)-1.
    flux_power = np.zeros(npix, dtype=np.double)
    # Calculate flux power for each spectrum in turn
    flux_power_perspectra = powerspectrum(dflux, axis=1)
    #Take the mean
    flux_power = np.mean(flux_power_perspectra, axis=0)
    assert np.shape(flux_power) == (npix//2+1,)
    return flux_power
