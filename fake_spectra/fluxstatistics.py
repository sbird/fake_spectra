# -*- coding: utf-8 -*-
"""Module to compute flux statistics from spectra:
the power spectrum, the pdf and to normalise to a mean tau.
Useful for lyman alpha forest work."""

import math
import os
import numpy as np
import scipy.fft
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor

# You need `nbodykit` only if you want to compute the 3D power spectrum with `flux_power_3d`
try :
    from nbodykit.lab import FFTPower
    from nbodykit.source.catalog import ArrayCatalog
    from nbodykit.source.mesh import ArrayMesh
    from nbodykit import setup_logging
    from nbodykit import CurrentMPIComm
except ImportError:
    pass

def obs_mean_tau(redshift):
    """The mean flux from 0711.1862: effective optical depth is (0.0023±0.0007) (1+z)^(3.65±0.21)
    Todo: check for updated values."""
    return 0.0023*(1.0+redshift)**3.65

#Batches smaller than this are not worth threading the FFT over.
_FFT_MINTHREAD = 1 << 20

#Threads for mean_flux, created on first use.
_MF_POOL = None
#Chunks smaller than this are not worth handing to another thread.
_MF_MINCHUNK = 250000

def _mean_flux_sums(tau, scale, out):
    """Partial sums of exp(-scale*tau) and tau*exp(-scale*tau) over one chunk of tau.
    The numpy ufuncs release the GIL, so chunks are summed in parallel."""
    flux = np.multiply(tau, -scale, out=out)
    np.exp(flux, out=flux)
    mean_flux = np.sum(flux)
    np.multiply(flux, tau, out=flux)
    return mean_flux, np.sum(flux)

def mean_flux(tau, mean_flux_desired, tol = 1e-5, nthreads=None):
    """Scale the optical depths by a constant value until we get the observed mean flux.
    ie, we want F_obs = bar{F} = < e^-tau >
    Solves iteratively using Newton-Raphson.
    This is safe because the exponential function is so well-behaved.
    Arguments:
        tau - optical depths to scale
        mean_flux_desired - mean flux desired
        tol - tolerance within which to hit mean flux
        nthreads - threads to use for the sums (default: all available cores)
    returns:
        scaling factor for tau"""
    global _MF_POOL
    tau = np.ravel(np.asarray(tau, dtype=np.float64))
    nbins = np.size(tau)
    if nbins == 0:
        return 0
    if nthreads is None:
        nthreads = len(os.sched_getaffinity(0))
    nchunk = max(1, min(nthreads, nbins // _MF_MINCHUNK))
    bounds = np.linspace(0, nbins, nchunk+1).astype(int)
    chunks = [tau[bounds[i]:bounds[i+1]] for i in range(nchunk)]
    #Scratch space, allocated once and reused by every iteration.
    scratch = [np.empty_like(cc) for cc in chunks]
    if nchunk > 1 and _MF_POOL is None:
        _MF_POOL = ThreadPoolExecutor(max_workers=nthreads)
    newscale = 1.
    while True:
        scale = newscale
        #Farm out all but the first chunk, then do that one here.
        futures = [_MF_POOL.submit(_mean_flux_sums, chunks[i], scale, scratch[i]) for i in range(1, nchunk)]
        sums = [_mean_flux_sums(chunks[0], scale, scratch[0])] + [ff.result() for ff in futures]
        flux = math.fsum([ss[0] for ss in sums])
        tau_flux = math.fsum([ss[1] for ss in sums])
        #Newton-Raphson
        newscale = scale + (flux - mean_flux_desired * nbins)/tau_flux
        #We don't want the absorption to change sign and become emission;
        #0 is too far.
        if newscale <= 0:
            newscale = 1e-10
        if abs(newscale - scale) <= tol * newscale:
            return newscale

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
    (nspec, npix) = np.shape(tau)
    mean_flux_power = np.zeros(npix//2+1, dtype=np.float64)
    #The k=0 mode of each sightline is the flux summed over pixels, which is
    #all we need to get the mean flux: no separate pass over tau required.
    kzero = np.empty(nspec, dtype=np.float64)
    # compute in batches, purely for computational efficiency
    for i in range(10):
        start = i*nspec//10
        end = min((i+1)*nspec//10, nspec)
        flux = np.exp(-scale*tau[start:end])
        # Calculate flux power for each spectrum in turn.
        # scipy's fft threads over the transforms, np.fft does not. The threads
        # spin between batches, so only ask for them when there is enough work.
        workers = -1 if flux.size >= _FFT_MINTHREAD else 1
        rfftd = scipy.fft.rfft(flux, axis=1, workers=workers, overwrite_x=True)
        kzero[start:end] = rfftd[:, 0].real
        mean_flux_power += np.sum(np.abs(rfftd)**2, axis=0)
    if mean_flux_desired is None:
        mean_flux_desired = np.sum(kzero)/(nspec*npix)
    #We want the power of d_F = F/mean(F) - 1. The FFT is linear, so dividing
    #by the mean flux just rescales every mode and subtracting one shifts k=0
    #alone: both can be applied to the summed power. The npix**2 normalises
    #the FFT so it is independent of input size, and vmax converts the units.
    mean_flux_power *= vmax/(npix**2 * nspec * mean_flux_desired**2)
    mean_flux_power[0] = vmax*np.sum((kzero/mean_flux_desired - npix)**2)/(npix**2 * nspec)
    mean_flux_power = mean_flux_power.astype(tau.dtype)
    assert np.shape(mean_flux_power) == (npix//2+1,)
    kf = _flux_power_bins(vmax, npix)
    #Divide out the window function
    if window and spec_res > 0:
        mean_flux_power /= _window_function(kf, R=spec_res, dv=vmax/npix)**2
    return kf,mean_flux_power

def _3d_powerspectrum(dflux_mesh, boxsize, los, dk=None, Nmu=10):
    """Compute the 3D power spectrum of the input using nbodykit
    Parameters:
    dfux_mesh - 3D array of flux variations, type is `mesh` in `nbodykit`
    boxsize - size of the box in units of interest (eg, comoving cMpc/h),
                the units of the 3d power spectrum, i.e. P(k,mu), will be in these units
    los - line of sight direction, i.e. [0,0,1] for z-axis
    dk - bin width in k
    Nmu - number of mu bins
    Returns:
    power - a dictionary with the p(k,mu) and the k and mu bins, keys:['power','k','mu']
    """
    power = FFTPower(dflux_mesh, BoxSize=boxsize,
                     mode='2d', los= los, dk=dk,
                     Nmu=Nmu)
    return power.power

def flux_power_3d(comm_nbodykit, tau, boxsize, mean_flux_desired=None, dk=None, Nmu=10, quiet=True):
    """Get the power spectrum of (variations in) the flux in 3D which is binned in (k,mu).
        This is: P_3D(k) = <d_F d_F>
                 d_F = e^-tau / mean(e^-tau) - 1
                 Then we bin P_3D(k) in k and mu.
        If mean_flux_desired is set, the spectral optical depths will be rescaled
        to match the desired mean flux.
        We compute the power spectrum along each sightline and then average the result.
        Arguments:
        comm_nbodykit: MPI communicator for nbodykit, I prefer to have one communicator for each process, i.e.
                        turning off parallel processing in nbodykit cause it is already fast enough.
                        You can set it as None if parallelism is not a concern to you.
            tau - optical depths. Shape is (NumLos, npix)
            mean_flux_desired - Mean flux to rescale to.
        boxsize - size of the box in units of interest (eg, comoving cMpc/h),
                the units of the 3d power spectrum, i.e. P(k,mu), will be in these units
        Returns:
            k, mu - the k and mu bins of the power spectrum
            flux_power - flux power spectrum in `boxsize` units
            Note: The first row corresponds to k=0, so you can remove it later
    """
    with CurrentMPIComm.enter(comm_nbodykit):
        scale = 1.
        if mean_flux_desired is not None:
            scale = mean_flux(tau, mean_flux_desired)
            tau *= scale
            print(f"rescaled: {scale}, mean_flux_desired = {mean_flux_desired}")
        else:
            mean_flux_desired = np.mean(np.exp(-tau))
        (nspec, npix) = np.shape(tau)
        nt = np.sqrt(nspec/3).astype(int)
        x, y, z = np.meshgrid(np.arange(nt), np.arange(nt), np.arange(npix), indexing='ij')
        x = x*boxsize/nt
        y = y*boxsize/nt
        z = z*boxsize/npix
        print(f'original dimenstions {(nt, nt, npix)}', flush=True)
        coords = np.vstack((x.ravel(), y.ravel(), z.ravel())).T
        for i in range(3):
            end = min((i+1)*nspec//3, nspec)
            # Turn on nbodkit's loging
            setup_logging('debug')
            # No interpoaltion is needed if the data is already on a uniform cube
            if npix == nt:
                mesh = ArrayMesh((np.exp(-tau[i*nspec//3:end])/mean_flux_desired -1 ).reshape((nt, nt, npix)), BoxSize=boxsize)
                mesh = mesh.compute(Nmesh=(nt,nt,nt))
            # Otherwise, do TSC interpoaltion to match the transverse resolution
            else:
                print(f'Interpolating the spectra along the perp direction | {datetime.now()}', flush=True)
                cat = ArrayCatalog({'Position': coords, 'df': np.exp(-tau[i*nspec//3:end].ravel()) / mean_flux_desired - 1})
                mesh = cat.to_mesh(Nmesh=[nt, nt, nt], value='df', BoxSize=boxsize, resampler='tsc', compensated=True, interlaced=True)

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
