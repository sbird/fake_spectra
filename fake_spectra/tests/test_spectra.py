# -*- coding: utf-8 -*-
"""
    Tests for the spectrum module.

Methods in spectra.py:

def __init__(self,num, base,cofm, axis, res=1., cdir=None, savefile="spectra.hdf5", savedir=None, reload_file = False, spec_res = 8):

#IO.
#Not tested.
def save_file(self):
def _save_file(self, f):
def _load_all_multihash(self,array, array_name):
def _save_multihash(self,save_array, grp):
def _really_load_array(self, key, array, array_name):
def load_savefile(self,savefile=None):

#Spectrum and array generation: these mostly either load from the snapshot or are done in C.
#Not tested here.
def _interpolate_single_file(self,fn, elem, ion, ll, get_tau):
def _read_particle_data(self,fn, elem, ion, get_tau):
def _filter_particles(self, elem_den, pos, velocity, den):
def _get_elem_den(self, elem, ion, den, temp, data, ind, ind2, star):
def _do_interpolation_work(self,pos, vel, elem_den, temp, hh, amumass, line, get_tau):
def particles_near_lines(self, pos, hh,axis=None, cofm=None):
def _vel_single_file(self,fn, elem, ion):
def _temp_single_file(self,fn, elem, ion):
def compute_spectra(self,elem, ion, ll, get_tau):
def get_mass_frac(self,elem, data, ind):
def replace_not_DLA(self, thresh=10**20.3):
def get_cofm(self, num = None):
def filter_DLA(self, col_den, thresh=10**20.3):
#This is unused and may not work.
def get_particle_number(self, elem, ion, res=8):

#Accessor methods.
#Not tested. Simple.
def get_tau(self, elem, ion,line, number = -1, force_recompute=False, noise=True):
def get_velocity(self, elem, ion):
def get_temp(self, elem, ion):
def get_col_density(self, elem, ion, force_recompute=False):
def get_metallicity(self):
def get_ion_metallicity(self, species,ion):
#Does non-trivial unit manipulation, should be tested eventually.
def get_density(self, elem, ion, force_recompute=False):
def equivalent_width(self, elem, ion, line):

#Connected to halo finding routines
#Not yet tested
def get_spectra_proj_pos(self, cofm=None):
def min_halo_mass(self, minpart = 400):
def load_halo(self):
def assign_to_halo(self, zpos, halo_radii, halo_cofm):
def find_nearest_halo(self):
def find_nearby_halos(self):
def get_contiguous_regions(self, elem="H", ion = 1, thresh = 2e20, relthresh = 1e-3):
def combine_regions(condition, mindist=0):
def contiguous_regions(condition):
def mass_hist(self, dm=0.1):
def virial_vel(self, halos=None, subhalo=False):

#Does nothing: just there to be overridden
def get_filt(self, elem, ion, thresh = 1):

#Transform the spectra for analysis.
#To test.
def add_noise(self, snr, tau, seed):

#Statistics.
#To test.
def _rho_abs(self, thresh=10**20.3, upthresh=10**40, elem = "H", ion = 1):
def rho_DLA(self, thresh=10**20.3):
def omega_abs(self, thresh=10**20.3, upthresh=10**40, elem = "H", ion = 1):
def omega_abs_cddf(self, thresh=10**20.3, upthresh=10**40, elem = "H", ion = 1):
def line_density(self, thresh=10**20.3, upthresh=10**40, elem = "H", ion = 1):
def line_density_eq_w(self, thresh=0.4, elem = "H", ion = 1, line=1216):
def column_density_function(self,elem = "H", ion = 1, dlogN=0.2, minN=13, maxN=23.):
def eq_width_hist(self, elem, ion, line, dv=0.05, eq_cut = 0.02):

"""

import numpy as np

from fake_spectra import spectra as ss
from fake_spectra import unitsystem
from fake_spectra import spec_utils
from fake_spectra import voigtfit
from fake_spectra import halocat

#def setup():
    #"""Load the fake data section and module to be used by these tests"""

def testRhoCrit():
    """Critical density at z=0"""
    units = unitsystem.UnitSystem()
    assert units.rho_crit(0.7) == 9.204285430050004e-30
    assert units.rho_crit(1.0) == 1.8784255979693885e-29

def testAbsDist():
    """Check absorption distance computation"""
    units = unitsystem.UnitSystem()
    assert units.absorption_distance(25000, 3) == 0.13377926628219666
    assert units.absorption_distance(25000, 2) == 0.07525083728373562
    assert units.absorption_distance(25000, 3) / units.absorption_distance(12500, 3) == 2.

def testRolledSpectra():
    """Check that we can correctly rotate a spectrum so the maximum is in the middle"""
    tau = np.zeros((2,50))
    tau[0,0] = 1
    tau[1,0] = 1
    tau[1,-1] = 2
    (roll, tau_new) = spec_utils.get_rolled_spectra(tau)
    assert np.all(roll == np.array([25,-24]))
    assert tau_new[0,25] == 1
    assert tau_new[1,25] == 2
    assert tau_new[1,26] == 1
    assert np.sum(np.abs(tau_new)) == 4

def testrescorr():
    """Check that convolving a spectrum with a Gaussian works"""
    tau = np.zeros((2,50))
    tau[0,25] = 2
    tau[1,23] = 3
    tau2 = spec_utils.res_corr(tau, 2, 8)
    #Check flux conserved
    assert np.abs(np.sum(tau2[0,:])/ np.sum(tau[0,:]) -1) < 1e-6
    assert np.abs(np.sum(tau2[1,:])/ np.sum(tau[1,:]) -1) < 1e-6
    #Check expanded by expected amount
    for i in (0,1):
        assert np.size(np.where(tau2[i,:]> 0)) == 15

# def testNan():
#     """Test nan"""
#     spec = ss.Spectra(3,'/home/spb/data/Cosmo/Cosmo0_V6/L25n512/output', cofm=np.array([[ 10724.84151495,   4444.02494373,  10534.57817268]]), axis=np.array([1,]), savefile="testfile.hdf5", reload_file=True)
#     tau = spec.get_tau("H", 1, 1215)
#     assert not np.any(np.isnan(tau))

def test_voigtfit():
    """Simple tests that the Voigt fitter is working, using 10 test CIV spectra."""
    import os.path

    fn = os.path.join(os.path.dirname(__file__), "example_civ_tau.npz")
    taus = np.load(fn)["arr_0"]
    for tau in taus:
        assert np.shape(tau) == (473,)
        prof = voigtfit.Profiles(tau,5.0103430332365999,elem="C",ion=4,line=1548)
        prof.do_fit()
        (ll, tfit) = prof.get_fitted_profile()
        #Check the fit is reasonable
        assert np.sum((tfit - tau)**2/(tau+0.5)**2)/np.size(tfit) < 0.05

def test_hcd_isolated():
    """Check that a single isolated HCD is fitted out and its core masked."""
    nbins = 2048
    dvbin = 5.
    prof = voigtfit.HCDProfiles(nbins, dvbin)
    amplitude = 5e6
    tau = np.roll(prof.profile(prof.btherm, amplitude), 300)
    (newtau, mask) = prof.do_hcd_fit(tau)
    #The wings should be subtracted to nothing: the residual is 5e-3
    assert np.max(newtau) < 1.
    #The core, where the profile is saturated, should be masked
    assert np.all(mask == (tau > 1))
    #The peak of the profile is masked, the edges of the spectrum are not
    assert mask[(nbins//2 + 300) % nbins]
    assert not mask.all()

def test_hcd_amplitude():
    """Check that the fitted amplitude of an isolated HCD is recovered."""
    prof = voigtfit.HCDProfiles(2048, 5.)
    for amplitude in (1e6, 5e6, 1e7, 1e8):
        tau = prof.profile(prof.btherm, amplitude)
        (_, peak_index, fitted) = prof.iterate_new_spectrum(tau)
        assert peak_index == 1024
        assert np.abs(fitted / amplitude - 1) < 0.01

def _forest(nbins, seed=7):
    """A smooth, noiseless forest-like optical depth field with mean tau ~0.4."""
    rng = np.random.default_rng(seed)
    kk = np.fft.rfftfreq(nbins)
    gauss = np.fft.irfft(np.fft.rfft(rng.normal(size=nbins))*np.exp(-0.5*(kk/0.01)**2), nbins)
    return np.exp(gauss/np.std(gauss) - 1.)

def test_hcd_forest():
    """Check that absorption blended with the HCD does not bias the fitted amplitude.
    The envelope estimator is used precisely because a least squares fit is biased
    high here by a factor of a few: extra absorption is one-sided, so the only way
    to account for it in a fit is to deepen the profile."""
    nbins = 2048
    prof = voigtfit.HCDProfiles(nbins, 5.)
    amplitude = 8e6
    hcd = prof.profile(prof.btherm, amplitude)
    #A uniform absorbing floor, and a smoothly varying forest
    for extra in (0.1, 0.3, 1.0, _forest(nbins)):
        fitted = prof.iterate_new_spectrum(hcd + extra)[2]
        assert np.abs(fitted / amplitude - 1) < 0.05
    #End to end: the forest outside the mask must survive the subtraction intact
    forest = _forest(nbins)
    (newtau, mask) = prof.do_hcd_fit(hcd + forest)
    #The mask should cover the core and not much more
    assert (hcd > 1).sum() <= mask.sum() < 1.1 * (hcd > 1).sum()
    #Nothing outside the mask should be over-subtracted down to zero optical depth
    assert not np.any((newtau[~mask] == 0) & (forest[~mask] > 0.1))
    assert np.abs(np.mean(np.exp(-newtau[~mask])) / np.mean(np.exp(-forest[~mask])) - 1) < 0.01

def test_hcd_saturated():
    """Check the fallback for an absorber saturated over the whole spectrum,
    so that there is no unsaturated shoulder to fit to."""
    prof = voigtfit.HCDProfiles(64, 5.)
    amplitude = 1e8
    tau = prof.profile(prof.btherm, amplitude)
    #Even the least absorbed pixel in the spectrum is saturated
    assert np.min(tau) > 100
    assert np.abs(prof.iterate_new_spectrum(tau)[2] / amplitude - 1) < 0.05

def test_hcd_blended():
    """Check that two blended HCDs are both found, and that a spectrum
    without an HCD is left alone."""
    nbins = 2048
    prof = voigtfit.HCDProfiles(nbins, 5.)
    rng = np.random.default_rng(23)
    noise = rng.exponential(0.3, nbins)
    tau = noise + np.roll(prof.profile(prof.btherm, 3e6), -600) + prof.profile(prof.btherm, 8e6)
    (newtau, mask) = prof.do_hcd_fit(tau)
    #Both cores are masked
    assert mask[nbins//2]
    assert mask[(nbins//2 - 600) % nbins]
    #Both damping wings are gone: what is left is the noise
    assert np.max(newtau) < 10 * np.max(noise)
    #A spectrum with no HCD in it should be untouched
    (newtau, mask) = prof.do_hcd_fit(noise)
    assert not np.any(mask)
    assert np.all(newtau == noise)
