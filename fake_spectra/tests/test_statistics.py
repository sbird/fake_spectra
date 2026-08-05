"""Module to test the flux statistics computation"""

import math
import numpy as np

from fake_spectra import fluxstatistics as stat

def testMeanFlux():
    """Test that we scale for the mean flux correctly"""
    tol = 1e-4
    nn = np.arange(1,101)
    #Use log x so that the mean flux is x^(-n) and stat.mean_flux returns n.
    tau = np.log(nn)
    mf = np.mean(np.exp(-tau))
    assert abs(stat.mean_flux(tau, mf,tol) - 1) < tol
    mf2 = np.mean(nn**(-2.))
    assert abs(stat.mean_flux(tau, mf2,tol) - 2) < tol
    mf3 = np.mean(nn**(-0.5))
    assert abs(stat.mean_flux(tau, mf3,tol) - 0.5) < tol

def testCalcPdf():
    """Test that we calculate the pdf of the flux correctly"""
    nn = np.arange(1,101,dtype=np.double)
    tau = np.log(nn)
    (bins,hist) = stat.flux_pdf(tau, 20)
    print(bins)
    assert bins[0] == 0.+1/40.
    assert bins[-1] == 1.-1./40.
    assert np.min(hist) == 0.
    assert np.max(hist) > 1.
    print(hist)
    expected = np.array([ 16. ,   2.2,   0.6,   0.2,   0.2,   0.2,   0.2,   0. ,   0. , 0. ,   0.2,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. , 0. ,   0.2])
    assert np.abs(np.sum(expected) - np.sum(hist)) < 1e-3
    #One of the points moves around from bin 1 to bin 2, depending on roundoff.
    assert np.all(np.abs(hist[3:] - expected[3:]) < 1e-3)
    assert np.abs(hist[0] - expected[0]) < 1e-2

def _sine_tau(npix, freq, amp, phase=0.):
    """Optical depths whose flux is exactly 1 + amp*sin(2 pi freq x), so that
    the mean flux is one and the flux variation d_F is a pure sine wave of
    known amplitude. freq must be a whole number of cycles across the box."""
    xx = np.arange(npix)/npix
    return -np.log1p(amp*np.sin(2*math.pi*(freq*xx + phase)))

def _parseval_sum(power, npix):
    """Sum a flux power spectrum over every DFT mode. The rfft drops the
    negative frequencies, which duplicate all but k=0 and, for an even number
    of pixels, the Nyquist mode."""
    total = 2*np.sum(power) - power[0]
    if npix % 2 == 0:
        total -= power[-1]
    return total

def testFluxPowerSine():
    """A flux variation which is a pure sine wave puts all the power in one
    bin, with an amplitude we know analytically."""
    amp = 0.4
    freq = 50
    #Check both even and odd binning
    for npix in (200, 201):
        taus = np.vstack([_sine_tau(npix, freq, amp),]*10)
        bins, power = stat.flux_power(taus, vmax=1., spec_res=0)
        assert np.shape(power) == (npix//2+1,)
        #A sine of amplitude a has power a^2/4 in its own bin and nothing anywhere else.
        assert np.argmax(power) == freq
        assert abs(power[freq] - amp**2/4) < 1e-12
        #The bins must line up with the power: the peak is at the frequency we
        #put in, in the 2 pi k / vmax convention flux_power returns.
        assert np.shape(bins) == np.shape(power)
        assert abs(bins[freq] - 2*math.pi*freq) < 1e-10
        assert bins[0] == 0.
        assert np.max(np.abs(np.delete(power, freq))) < 1e-20
        #Every sightline has the mean flux, so there is nothing in the k=0 mode.
        assert power[0] < 1e-20
        #vmax just sets the units
        _, power2 = stat.flux_power(taus, vmax=7., spec_res=0)
        assert np.all(np.abs(power2 - 7*power) < 1e-12)

def testFluxPowerParseval():
    """Check the normalisation: summed over modes the power is the variance
    of the flux variation."""
    rng = np.random.default_rng(23)
    for npix in (200, 201):
        taus = np.abs(rng.standard_normal((5, npix)))
        _, power = stat.flux_power(taus, vmax=1., spec_res=0)
        flux = np.exp(-taus)
        dflux = flux/np.mean(flux) - 1.
        assert abs(_parseval_sum(power, npix) - np.mean(dflux**2)) < 1e-14

def testFluxPowerMean():
    """Check that we take the mean of the power over sightlines. Each sine has
    the same mean flux, so they can be averaged without renormalising."""
    npix = 200
    tau1 = _sine_tau(npix, 50, 0.4)
    tau2 = _sine_tau(npix, 20, 0.25)
    _, p1 = stat.flux_power(np.vstack([tau1,]), vmax=1., spec_res=0)
    _, p2 = stat.flux_power(np.vstack([tau2,]), vmax=1., spec_res=0)
    _, both = stat.flux_power(np.vstack([tau1, tau2]), vmax=1., spec_res=0)
    assert np.all(np.abs(both - (p1+p2)/2) < 1e-14)
    #More sightlines than the batch loop has batches, with random phases
    rng = np.random.default_rng(5)
    taus = np.vstack([_sine_tau(npix, 50, 0.4, phase=pp) for pp in rng.random(25)])
    _, many = stat.flux_power(taus, vmax=1., spec_res=0)
    #The power does not care about the phase, so this is the single sightline answer
    assert np.all(np.abs(many - p1) < 1e-12)

def testFluxPowerWindow():
    """Check that the window function is divided out of the power"""
    npix = 200
    taus = np.vstack([_sine_tau(npix, 50, 0.4),]*10)
    bins, wind_power = stat.flux_power(taus, vmax=1., spec_res=0.01, window=True)
    _, power = stat.flux_power(taus, vmax=1., spec_res=0.01, window=False)
    wind = stat._window_function(bins, R=0.01, dv=1./npix)
    assert np.all(np.abs(wind_power*wind**2 - power) < 1e-12)

def testFluxPowerRescaled():
    """Check the power of optical depths rescaled to a desired mean flux."""
    rng = np.random.default_rng(24)
    mean_flux_desired = 0.3
    for npix in (200, 201):
        taus = np.abs(rng.standard_normal((5, npix)))
        _, power = stat.flux_power(taus, vmax=1., spec_res=0, mean_flux_desired=mean_flux_desired)
        scale = stat.mean_flux(taus, mean_flux_desired)
        flux = np.exp(-scale*taus)
        assert abs(np.mean(flux) - mean_flux_desired) < 1e-5*mean_flux_desired
        dflux = flux/mean_flux_desired - 1.
        assert abs(_parseval_sum(power, npix) - np.mean(dflux**2)) < 1e-12
        #Rescaling changed the answer
        _, unscaled = stat.flux_power(taus, vmax=1., spec_res=0)
        assert not np.allclose(power, unscaled)
