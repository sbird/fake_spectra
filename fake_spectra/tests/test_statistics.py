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
    nn = np.arange(1,101)
    tau = np.log(nn)
    (bins,hist) = stat.flux_pdf(tau, 20)
    print(bins)
    assert bins[0] == 0.+1/40.
    assert bins[-1] == 1.-1./40.
    assert np.min(hist) == 0.
    assert np.max(hist) > 1.
    assert np.all(np.abs(hist - np.array([ 16. ,   2.2,   0.6,   0.2,   0.2,   0.2,   0.2,   0. ,   0. , 0. ,   0.2,   0. ,   0. ,   0. ,   0. ,   0. ,   0. ,   0. , 0. ,   0.2])) < 1e-3)

def testPowerspectrum():
    """Test the core FFT routine"""
    xx = np.tile(np.arange(0,1,0.01)**2,(10,1))
#     xx = np.arange(0,1,0.01)**2
    assert np.shape(xx)[0] == 10
    fpk = stat._powerspectrum(np.exp(-xx), axis=1)
    #Check that Parseval's theorem is true, accounting for the negative frequency modes not included in the DFT.
    for ff in fpk:
        dpower = (np.sum(ff)+np.sum(ff[1:]))
        assert abs(dpower - np.sum(np.exp(-xx[0,:])**2)) < 0.1
    xx = np.linspace(1,51,200)
    inn = np.sin(2*math.pi*xx)
    #This will be a delta function centered at 50
    ff = stat._powerspectrum(inn)
    assert np.where(np.max(ff) == ff)[0][0] == 50

def testFluxPower():
    """Check that taking the mean of the flux power spectra is working correctly"""
    #Check both even and odd binning
    for bb in (200, 201):
        xx = np.linspace(0,51,bb)
        inn = np.sin(2*math.pi*xx)+1.5
        ff = stat._powerspectrum(np.exp(-inn)-1)
        #Construct some optical depths offset from each other
        taus = np.zeros((9,bb))
        taus = np.vstack([inn,]*10)
        bins, power = stat.flux_power(taus, 1.)
        power/=12.5569
        assert np.all(np.abs(power[1:] - ff[1:]) < 0.01*ff[1:])
        assert power[0] < 1e-20
