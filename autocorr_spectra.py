# -*- coding: utf-8 -*-
"""Module to compute the autocorrelation of quantities along spectra"""
import _autocorr_spectra_priv


def autocorr_spectra(slist, spos, pixsz, nbins=100):
    """
    Find the autocorrelation function from a list of spectra
    Spectra are assumed to be along the same axis.
    slist - list of quantity along spectra to autocorrelate. npix * nspectra
    spos -  positions of the spectra: 2x nspectra: (x, y).
    nbins - number of bins in output autocorrelation function
    pixsz - Size of a pixel in units of the spectra position.
    """
    (modes, auto) = _autocorr_spectra_priv.autocorr_spectra(slist, spos, pixsz, nbins)
    auto /= modes
    return auto

