# -*- coding: utf-8 -*-
"""Contains the power spectrum - specific functions for the spectrum analysis code."""

from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt

from . import spectra
from . import fluxstatistics as fstat

try:
    xrange(1)
except NameError:
    xrange = range


class CalcPowerspectrum(spectra.Spectra):
    """Class to calculate power spectrum and associated things.
    I assume we already have a spectra file written."""

    def __init__(self, num, base, savefile, **kwargs):
        spectra.Spectra.__init__(
            self, num, base, cofm=None, axis=None, savefile=savefile, reload_file=False, **kwargs)

    def calc_scaling(self, elem="H", ion=1, line=1215, mean_flux_desired=None):
        """Calculate the scaling factor of UVB in order to obtain mean_flux_desired"""
        if mean_flux_desired is None:
            return 1.

        tau = self.get_tau(elem, ion, line)
        return fstat.mean_flux(tau, mean_flux_desired)

    def calc_powerspectrum(self, elem="H", ion=1, line=1215, mean_flux_=None, mean_flux_desired=None, scale=1., window=True, kmin=0.001, kmax=0.1, N=1000):
        """Calculate power spectrum and do an interpolation in the range [kmin, kmax], using N points. Return the interpolated values.
        Can choose to input mean_flux_desired or scaling of the UVB."""
        xinterp = np.linspace(np.log10(kmin), np.log10(kmax), N)
        rst = self.get_flux_power_1D(
            elem, ion, line, mean_flux_, mean_flux_desired, scale, window)
        yinterp = np.interp(xinterp, np.log10(rst[0]), np.log10(
            rst[1]), left=np.nan, right=np.nan)
        return 10**xinterp, 10**yinterp
