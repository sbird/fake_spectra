"""
A simple Voigt profile fitter. Works by finding peaks and removing them iteratively.
All the peaks are then re-fit to the spectrum at once. Minimisation is done by scipy.
Algorithm suggested by Boris Leistedt, closely based on that in AUTOVP by Ben Oppenheimer & Romeel Dave.
"""

import math
import time
import multiprocessing
import numpy as np
import scipy.optimize as optimize
import scipy.special

from . import line_data

class Profiles(object):
    """Class to fit a spectrum input (without noise) with a series of profiles, Gaussian or Voigt.
    Inputs: tau - spectrum to fit.
            dvbin - Size of velocity bin in km/s.
            profile - either Voigt of Gaussian.
            elem, ion, line - Line to fit and ion to use.
    """
    def __init__(self, tau, dvbin, profile="Voigt", elem="H", ion=1, line=1215):
        self.dvbin = dvbin
        self.amplitudes = []
        self.means = []
        self.stddev = []
        if profile == "Gaussian":
            self.profile = self.gaussian_profile
        else:
            self.profile = self.voigt_profile
        lines = line_data.LineData()
        line = lines[(elem, ion)][line]
        #The /10^13 converts from Angstrom to km, so this is km/s
        self.voigt_fac = line.gamma_X*line.lambda_X/(4*math.pi)/1e13
        sigma_t = 6.652458558e-25 # Thompson cross-section in cm^2
        #lambda_x is in Angstrom, want cm
        self.sigma_a = math.sqrt(3.0*math.pi*sigma_t/8.0) * line.lambda_X/1e8  * line.fosc_X
        #Speed of light in km/s
        self.light = 2.99e5
        #Set things up
        self.tau = tau
        #In km/s (the units of dvbin)
        self.wavelengths = np.arange(0, np.size(tau))*self.dvbin
        #Wavelength difference from central region
        midpt = np.size(self.wavelengths)//2
        self.lambda_diff = (self.wavelengths - midpt*self.dvbin)

    def do_fit(self, tol=1e-4, signif=0.95, sattau=4):
        """Do the fit.
        tol - stop peak finding when the largest peak is tol * max(1,max(self.tau))
        As we iteratively remove peaks, small peaks are increasingly likely to just be fitting errors
        from the larger ones. Default value is 10^-4 because we cut off the tails of the Gaussian at 10^-5.
        signif - We will do a global re-fit of all the peaks at the end. If adding the extra peak
        improves the fit by less than this value, the remaining peaks are considered not significant.
        """
        fitted_tau = np.array(self.tau)
        #We do not care about very small peaks; this includes both peaks which are small in absolute value
        #and peaks which are small in a relative sense.
        realtol = tol * np.max([1, np.max(self.tau)])
        #Initialize mask using the whole spectrum for fitting the first peak
        mask = np.where(fitted_tau > -1)
        #Do the fit iteratively, stopping when the largest peak is likely unimportant.
        while mask[0].size != 0 and np.max(fitted_tau[mask]) > realtol:
            (fitted_tau, mean, amplitude, stddev) = self.iterate_new_spectrum(fitted_tau, mask=mask)
            #We only want to fit to regions that are not already saturated in the fit.
            self.means.append(mean)
            self.amplitudes.append(amplitude)
            self.stddev.append(stddev)
            mask = np.where(self.profile_multiple(self.stddev, self.means, self.amplitudes) < sattau)
        #Do a global re-fit of all peaks
        inputs = np.hstack([self.stddev[:2], self.means[:2], self.amplitudes[:2]])
        # tolerance for global re-fit should be larger than for single profiles -- set to 25% of largest peak
        # or 0.25, whichever is greater
        result = optimize.minimize(self.fun_min_multiple, inputs, tol=realtol/tol*0.25, method='Nelder-Mead')
        total = np.min([2, np.size(self.stddev)])
        prior_result = result.fun*2
        for maxpk in range(3, np.size(self.stddev)+1):
            inputs = np.hstack([self.stddev[:maxpk], self.means[:maxpk], self.amplitudes[:maxpk]])
            new_result = optimize.minimize(self.fun_min_multiple, inputs, tol=realtol/tol*0.25, method='Nelder-Mead')
            #If the addition of the previous two peaks does not substantially improve the fit, stop.
            if new_result.fun/prior_result > signif:
                total = maxpk-1
                break
            #If it did improve the fit, but didn't converge, stop anyway.
            if not new_result.success:
                total = maxpk
                result = new_result
                break
            #Otherwise, add the extra peak and continue
            total = maxpk
            prior_result = result.fun
            result = new_result
        assert np.size(result.x) == total*3
        self.stddev_new = np.abs(result.x[0:total])
        self.means_new = result.x[total:2*total] % np.max(self.wavelengths)
        self.amplitudes_new = result.x[2*total:]
#         refitted=self.profile_multiple(self.stddev_new, self.means, self.amplitudes)
#         assert np.all(np.abs((tol+refitted)/(tol+self.tau) -1.) < 0.5)

    def iterate_new_spectrum(self, f, mask=None):
        """Fit the largest peak and return a new spectrum with that peak removed."""
        #Find largest peak in unmasked region
        if mask is None:
            peak_index = np.argmax(f)
        else:
            peak_index = np.where(f == np.max(f[mask]))[0][0]
        mean = peak_index*self.dvbin
        amplitude = f[peak_index]
        #First roll the spectrum to avoid edge effects
        midpt = np.size(f)//2
        maxx = midpt - peak_index
        f_rolled = np.roll(f, maxx)
        #Do the fit for the width
        optargs = (amplitude, f_rolled)
        result = optimize.minimize_scalar(self.fun_min, bracket=(10., 120.), bounds=(1., 120.), method='bounded', args=optargs)
        fitted = self.profile(result.x, mean, amplitude)
        assert np.argmax(fitted) == peak_index
        newf = f-fitted
        assert np.max(newf) < np.max(f)
        return newf, mean, amplitude, result.x

    def fun_min(self, stddev, amplitude, tau):
        """Helper function to pass to scipy.optimise. Computes the differences
        between the profile and the input spectrum.
        As each spectrum should be localised, down-weight far-off points.
        Function is assumed to be already rotated so that the max value is in the middle."""
        midpt = np.size(tau)//2
        assert tau[midpt] == amplitude
        gauss = self.profile(stddev, midpt*self.dvbin, amplitude)
        assert np.argmax(gauss) == midpt
        #Try to avoid separated peaks messing with the fit.
        #We find a minima between the peak and each edge, and exclude points beyond it from the fit.
        ind = np.arange(np.size(tau)-2)+1
        mins = np.where(np.logical_and(tau[ind] <= tau[ind+1], tau[ind] <= tau[ind-1]))[0]
        #Check whether the edges are a minima.
        #This is more obvious to read than mucking about with %
        if tau[-1] < tau[-2] and tau[-1] < tau[0]:
            mins = np.append(mins, -1)
        if tau[0] < tau[-1] and tau[0] < tau[1]:
            mins = np.append(mins, 0)
        assert np.size(mins) > 0
        # only use minima that are low enough compared to peak (~0.08 different in flux)
        mins = mins[np.where(np.abs(amplitude - tau[mins]) > np.log(1./0.92))]
        minn = np.max(np.append(mins[np.where(mins < midpt)], 0))
        maxx = np.min(np.append(mins[np.where(mins > midpt)], np.size(tau)))
        assert minn < midpt and maxx > midpt
        diff = np.sum(((np.exp(-tau) - np.exp(-gauss)))[minn:maxx]**2)
        return diff

    def profile_multiple(self, stddevs, means, amplitudes):
        """Function to generate a profile consisting of multiple peaks, so we can re-fit one last time."""
        gauss = np.sum([self.profile(stddev, mean, amplitude) for (stddev, mean, amplitude) in zip(stddevs, means, amplitudes)], axis=0)
        return gauss

    def fun_min_multiple(self, inputs):
        """Helper function for scipy.optimise that trys to fit the whole (original) spectrum at once
        with the detected peaks."""
        third = np.size(inputs)//3
        stddevs = inputs[0:third]
        means = inputs[third:2*third]
        amplitudes = inputs[2*third:]
        gauss = self.profile_multiple(stddevs, means, amplitudes)
        return np.sum((np.exp(-self.tau) - np.exp(-gauss))**2)

    def voigt_profile(self, stddev, mean, amplitude):
        """Compute the Voigt profile, which is the real part of the
           Faddeeva (complex probability) function of the variable
           w = u + i a
           So that F(w) = H(a,u) + i J(a,u) = exp(-w^2) erfc(-iw)
           Arguments:
           T0 = delta_v/btherm (velocity difference from the mean)
           aa: voigt_fac/btherm
           voigt_fac = gamma * lambda/(4*pi)/1e5 (a property of the line)
           (note btherm is sqrt(2k T / M))
           stddev = btherm
           Normalise so that at x = mean, return amplitude.
        """
        #Normalise the input - it is easier than putting constraints on the solvers.
        stddev = np.abs(stddev)
        amplitude = np.abs(amplitude)
        peak_index = mean/self.dvbin
        if peak_index < 0:
            peak_index += np.size(self.wavelengths)-1
        peak_index = int(round(peak_index))
        peak_index = peak_index % np.size(self.wavelengths)
        aa = 1j*self.voigt_fac/stddev
        #This is the Fadeeva function
        T0 = self.lambda_diff / stddev + aa
        fadeeva = np.real(scipy.special.wofz(T0))
        norm = np.real(scipy.special.wofz(aa))
        #Note this normalisation means that the profile is normalised to unity at its peak,
        #rather than having unity integral.
        rolled = fadeeva/norm*amplitude
        profile = np.roll(rolled, peak_index-np.argmax(rolled))
        assert np.argmax(profile) == peak_index
        return profile

    def gaussian_profile(self, stddev, mean, amplitude):
        """A Gaussian profile shape"""
        midpt = np.size(self.wavelengths)//2
        gauss = np.exp(-(self.wavelengths-midpt)**2/stddev**2/2)*amplitude
        profile = np.roll(gauss, int(round(mean/self.dvbin)) - np.argmax(gauss))
        assert np.argmax(profile) == int(round(mean/self.dvbin))
        return profile

    def get_fitted_profile(self):
        """Helper function to return the fitted profile shape"""
        return (self.wavelengths, self.profile_multiple(self.stddev_new, self.means_new, self.amplitudes_new))

    def get_b_params(self):
        """Helper function to return the doppler b (in km/s): for the definition used b = sigma"""
        return np.abs(self.stddev_new)

    def get_positions(self):
        """Helper function to return the wavelength in km/s (offset from a 0 at the start of the box) of each absorber."""
        return np.array(self.means_new)

    def get_column_density(self, btherm, amp):
        """Helper function to return the column density in 1/cm^-2 given an amplitude and a dispersion."""
        #Our profile is normalised to have peak value of self.amplitudes
        #The usual normalisation is integral_R voigt  dv = 1.
        #Correction from this is: amplitudes * b * sqrt(pi) / W(i gamma/b)
        #So we have N sigma_a c = int_R tau dv
        # N = 1/(sigma_a c) amplitudes b sqrt(pi) / W(i gamma/b)
        # vnorm is a cross-section in cm^-2.
        vnorm = btherm/self.light * math.sqrt(math.pi) / self.sigma_a
        #This is the Fadeeva function normalisation.
        fnorm = np.real(scipy.special.wofz(1j*self.voigt_fac/btherm))
        #Find amplitude divided by Voigt profile, which is still dimensionless.
        colden = amp * vnorm / fnorm
        return colden

    def get_column_densities(self):
        """Helper function to return the column densities in 1/cm^-2 of each absorber."""
        return np.array([self.get_column_density(bb, A) for (bb, A) in zip(np.abs(self.stddev_new), self.amplitudes_new)])

    def get_systems(self, close):
        """Get the absorption components merged into systems. This basically runs friends of friends on each spectrum and
        merges the closest systems until the closest systems are more than close km/s apart.
        Only the column density and positions are defined for the merged systems; b parameter doesn't make sense."""
        pos = self.get_positions()
        sorter = np.argsort(pos)
        colden = self.get_column_densities()[sorter]
        self.stddev_new = self.stddev_new[sorter]
        pos = pos[sorter]
        #Get maximum and minimum extent of absorbers: this starts off as just their peaks.
        maxpos = np.array(pos)
        #Move the end of the largest object back before the beginning of the spectrum.
        maxpos[-1] -= np.max(self.wavelengths)
        minpos = np.roll(np.array(pos), -1)
        #Compute difference in position between adjacent objects
        #this should include distance between the final and initial items
        distances = minpos-maxpos
        assert np.all(distances > 0)
        #As long as some absorbers are close
        while np.min(distances) < close:
            minn = np.argmin(distances)
            #Add the column densities together, average positions
            if minn == np.size(distances)-1:
                minn = -1
                pos[minn+1] = (pos[minn]+pos[minn+1]+np.max(self.wavelengths))/2. % np.max(self.wavelengths)
            else:
                pos[minn+1] = (pos[minn]+pos[minn+1])/2.
            colden[minn+1] += colden[minn]
            colden = np.delete(colden, minn)
            pos = np.delete(pos, minn)
            maxpos = np.delete(maxpos, minn)
            minpos = np.delete(minpos, minn)
            if len(pos) < 2:
                break
            distances = minpos-maxpos
            assert np.all(distances > 0)
        return colden, pos


class _SingleProfileHelper(object):
    """Picklable helper class to Voigt fit a single profile and optionally print how long it took.
    Used because lambdas are not picklable and functools.partial is not picklable on python 2."""
    def __init__(self, dvbin, elem, ion, line, verbose=False, close=0.):
        self.dvbin = dvbin
        self.elem = elem
        self.ion = ion
        self.line = line
        self.verbose = verbose
        self.close = close

    def __call__(self, tau_t):
        """Call the fit"""
        stime = time.time()
        prof = Profiles(tau_t, self.dvbin, elem=self.elem, ion=self.ion, line=self.line)
        prof.do_fit()
        n_this, _ = prof.get_systems(self.close)
        ftime = time.time()
        if self.verbose:
            print("Fit: systems=", np.size(n_this), np.size(prof.get_b_params()), "N=", np.max(n_this))
            print("Fit took: ", ftime-stime, " s")
        return n_this, prof.get_b_params()

def get_voigt_fit_params(taus, dvbin, elem="H", ion=1, line=1215, verbose=False, close=0.):
    """Helper function to get the Voigt parameters, N_HI and b in a single call."""
    return get_voigt_systems(taus, dvbin, elem, ion, line, verbose, close, b_param=True)

def get_voigt_systems(taus, dvbin, elem="H", ion=1, line=1215, verbose=False, close=0., b_param=False):
    """Helper function to get the Voigt parameters, N_HI and (optionally) b in a single call."""
    start = time.time()
    #Set up multiprocessing pool: lambdas are not picklable, so not using them.
    #functools.partial not picklable on python 2.
    helper = _SingleProfileHelper(dvbin, elem, ion, line, verbose, close)
    pool = multiprocessing.Pool(None)
    results, b_results = zip(*pool.map(helper, taus))
    n_vals = np.concatenate(results)
    end = time.time()
    print("Total fit took: ", end-start, " s")
    if b_param:
        b_params = np.concatenate(b_results)
        return n_vals, b_params
    return n_vals

def _power_fit(ln, lb0, gamm1):
    """A power law fit given some parameters"""
    return lb0 + gamm1 * (ln - 13.6)

def _opt_power_fit(inputs, lb_params, ln_vals):
    """Fit log b = log b_0  + (Gamma-1) * (log NHI - log NHI,0)
    for log b_0 and Gamma-1. log NHI,0 = 13.6"""
    lb0 = inputs[0]
    gamm1 = inputs[1]
    return np.sum(abs(lb_params - _power_fit(ln_vals, lb0, gamm1)))

def get_b_param_dist(taus, dvbin, elem="H", ion=1, line=1215):
    """Get the power law betweeen the 'minimum' b parameter and column density,
    following Rudie 2012 and Schaye 1999."""
    n_vals, b_params = get_voigt_systems(taus, dvbin, elem=elem, ion=ion, line=line, verbose=False, close=-1., b_param=True)
    used = np.where((b_params > 8)*(b_params < 100)*(n_vals < 10**(14.5))*(n_vals > 10**(12.5)))
    sort = np.argsort(np.log10(n_vals[used]))
    ln_vals = np.log10(n_vals[used])[sort]
    lb_params = np.log10(b_params[used])[sort]

    ###################### sigma rejection ###########################
    # will add each bin to final_retained as we iterate through
    # this works because the values are sorted
    final_retained = np.empty(0, dtype=bool)

    for i in range(8):
        bin_edges = np.where((ln_vals >= 12.5+i*0.25)*(ln_vals < 12.5+(i+1)*0.25))
        # current bin
        bin_range = ((ln_vals >= 12.5+i*0.25)*(ln_vals < 12.5+(i+1)*0.25)*(10**lb_params < 40))
        n_removed = 1
        while n_removed > 0:
            dev = np.std(lb_params[bin_range]) # stddev in the bin, for this iteration
            mean = np.mean(lb_params[bin_range]) # mean in the bin, for this iteration
            # keep includes all points, but we will only add those in the bin to final_retained
            keep = (abs(lb_params - mean) < 2*dev)
            n_removed = np.sum(bin_range) - np.sum(bin_range*keep)
            # only points which have true for all test so far will still be true
            bin_range = bin_range*keep

        # switch previously omitted points which are above the mean to true, including points over 40
        bin_range[(bin_range is False)*(lb_params > mean)] = True
        final_retained = np.append(final_retained, bin_range[bin_edges[0][0]:bin_edges[0][-1]+1])

    ln_vals = ln_vals[final_retained]
    lb_params = lb_params[final_retained]
    #################################################################

    # The fit proceeds iteratively.
    result = np.array([np.log10(18), 0.17])
    # First fit log b = log b_0 + (Gamma -1) * (log NHI - log NHI,0)
    # for Gamma-1 and b_0 with NHI,0 = 13.6.
    newresult = optimize.minimize(_opt_power_fit, result, args=(lb_params, ln_vals), tol=1e-4).x
    n_discarded = n_vals.size - used[0].size
    while n_discarded > 0:
        result = newresult
        # Points with b > 1 sigma above the fit are discarded and the fit is made again
        # until convergence (i.e. no more points are removed in a step).
        residuals = lb_params - _power_fit(ln_vals, result[0], result[1])
        sigma_fit = np.sum(abs(np.mean(residuals) - residuals))/residuals.size
        #remove points more than 1 sigma *above* fit. Note absence of np.abs!
        used = np.where(residuals <= sigma_fit)
        n_discarded = ln_vals.size - used[0].size
        lb_params = lb_params[used]
        ln_vals = ln_vals[used]
        #Make fit
        newresult = optimize.minimize(_opt_power_fit, result, args=(lb_params, ln_vals), tol=1e-4).x
    # Now all points > 1 sigma *below* the fit are removed and the fit is made again.
    residuals = lb_params - _power_fit(ln_vals, result[0], result[1])
    sigma_fit = np.sum(abs(np.mean(residuals) - residuals))/residuals.size
    #remove points more than 1 sigma *below* fit. Note absence of np.abs!
    used = np.where(residuals > -1*sigma_fit)
    lb_params = lb_params[used]
    ln_vals = ln_vals[used]
    #Make fit
    newresult = optimize.minimize(_opt_power_fit, result, args=(lb_params, ln_vals), tol=1e-4).x
    #Return the new b_min and Gamma.
    return (10**newresult[0], newresult[1])
