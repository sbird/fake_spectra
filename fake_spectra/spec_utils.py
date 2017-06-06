"""Some utility functions for the spectra."""
import numpy as np
from scipy.ndimage.filters import gaussian_filter1d

def res_corr(flux, dvbin, fwhm=8):
    """
        Real spectrographs have finite spectral resolution.
        Correct for this by smoothing the spectrum (the flux) by convolving with a Gaussian.
        The input spectrum is assumed to have infinite resolution, since we have used a spline
        to interpolate it first and/or we are converged.
        Strictly speaking we should rebin the spectrum after to have the same resolution
        as the observed pixel, but as long as the pixels are smaller than the FWHM of the
        spectrograph (which is the case as long as the observer is smart) we will be fine.
        args:
            flux - The input flux spectra
            dvbin - the width in km/s for the input flux
            fwhm - FWHM of the spectrograph in km/s
    """
    # Convert FWHM input to internal units
    res = fwhm/dvbin
    #FWHM of a Gaussian is 2 \sqrt(2 ln 2) sigma
    sigma = res/(2*np.sqrt(2*np.log(2)))
    #Do filter in wrapping mode to avoid edge effects
    oflux = gaussian_filter1d(flux, sigma, axis=-1, mode='wrap')
    return oflux

def get_rolled_spectra(tau):
    """
    Cycle the array tau so that the peak is at the middle.
    Returns (roll - the index the array was rolled by, tau_out - the rolled array)
    """
    (roll, tau_out) = zip(*[_roll_one_spectra(tau_l) for tau_l in tau])
    assert np.all(np.shape(roll) == np.shape(tau[:,0]))
    assert np.all(np.shape(tau_out) == np.shape(tau))
    return (np.array(roll), np.array(tau_out))

def _roll_one_spectra(tau_l):
    """Roll a single spectrum so the peak is in the middle."""
    max_t = np.max(tau_l)
    ind_m = np.where(tau_l == max_t)[0][0]
    tau_out = np.roll(tau_l, int(np.size(tau_l)/2)- ind_m)
    roll = int(np.size(tau_l)/2) - ind_m
    return roll, tau_out
