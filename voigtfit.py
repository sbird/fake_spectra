"""A simple Voigt profile fitter. Works by finding peaks and removing them iteratively.
Algorithm suggested by Boris Leistedt."""

import math
import numpy as np
import scipy.optimize as optimize
import scipy.special
import line_data

class Profiles(object):
    """Class to fit a spectrum input (without noise) with a series of profiles, Gaussian or Voigt.
    Inputs: tau - spectrum to fit.
            dvbin - Size of velocity bin in km/s.
            profile - either Voigt of Gaussian.
            elem, ion, line - Line to fit and ion to use.
            minweight - when fitting single profiles, the least squares is weighted by a gaussian.
                This is the std. dev, so points further than 2*distweight are not really considered.
                In km/s"""
    def __init__(self, tau, dvbin, profile="Voigt", elem="H", ion=1, line=1215, minweight=50.):
        self.dvbin = dvbin
        self.distweight = minweight
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
        self.wavelengths = np.arange(0,np.size(tau))*self.dvbin

    def do_fit(self, tol=1e-4, signif = 0.9):
        """Do the fit.
        tol - stop peak finding when the largest peak is tol * max(self.tau)
        As we iteratively remove peaks, small peaks are increasingly likely to just be fitting errors
        from the larger ones.
        signif - We will do a global re-fit of all the peaks at the end. If adding the extra peak
        improves the fit by less than this value, the remaining peaks are considered not significant.
        """
        fitted_tau = np.array(self.tau)
        #Do the fit iteratively, stopping when the largest peak is likely unimportant.
        while np.max(fitted_tau) > tol*np.max(self.tau):
            (fitted_tau, mean, amplitude, stddev)=self.iterate_new_spectrum(fitted_tau)
            self.means.append(mean)
            self.amplitudes.append(amplitude)
            self.stddev.append(stddev)
        #Do a global re-fit of all peaks
        total = np.size(self.stddev)
        inputs = np.hstack([self.stddev[:2], self.means[:2], self.amplitudes[:2]])
        result = optimize.minimize(self.fun_min_multiple,inputs, tol=tol*np.max(self.tau), method='Powell')
        for maxpk in range(3,total):
            inputs = np.hstack([self.stddev[:maxpk], self.means[:maxpk], self.amplitudes[:maxpk]])
            new_result = optimize.minimize(self.fun_min_multiple,inputs, tol=tol*np.max(self.tau), method='Powell')
            if not new_result.success:
                raise RuntimeError(new_result.message)
            #If the extra peak does not substantially improve the fit, stop.
            if new_result.fun > result.fun*signif:
                total = maxpk-1
                break
            result = new_result
        assert np.size(result.x) == total*3
        self.stddev_new = result.x[0:total]
        self.means_new = result.x[total:2*total]
        self.amplitudes_new = result.x[2*total:]
#         refitted=self.profile_multiple(self.stddev_new, self.means, self.amplitudes)
#         assert np.all(np.abs((tol+refitted)/(tol+self.tau) -1.) < 0.5)

    def iterate_new_spectrum(self, f):
        """Fit the largest peak and return a new spectrum with that peak removed."""
        #Find largest peak
        peak_index = np.argmax(f)
        mean = peak_index*self.dvbin
        amplitude = f[peak_index]
        #First roll the spectrum to avoid edge effects
        midpt = int(round(np.size(f)/2))
        maxx = midpt - peak_index
        f_rolled = np.roll(f,maxx)
        #Do the fit for the width
        optargs = (amplitude, f_rolled)
        result = optimize.minimize_scalar(self.fun_min,bracket=(10.,100.),bounds=(1.,100.),method='bounded',args=optargs)
        #The documentation says this is how you test success, but I guess it lies!
        #if not result.success:
        #    raise RuntimeError(result.mesg)
        return (f-self.profile(result.x,mean,amplitude), mean, amplitude, result.x)

    def fun_min(self, stddev, amplitude, tau):
        """Helper function to pass to scipy.optimise. Computes the differences
        between the profile and the input spectrum.
        As each spectrum should be localised, down-weight far-off points.
        Function is assumed to be already rotated so that the max value is in the middle."""
        midpt = int(round(np.size(tau)/2))
        assert np.argmax(tau) == midpt
        gauss = self.profile(stddev,midpt*self.dvbin, amplitude)
        assert np.argmax(gauss) == midpt
        #Try to avoid separated peaks messing with the fit.
        #We find a minima between the peak and each edge, and exclude points beyond it from the fit.
        ind = np.arange(np.size(tau)-2)+1
        mins = np.where(np.logical_and(tau[ind] <= tau[ind+1], tau[ind] <= tau[ind-1]))[0]
        #Check whether the edges are a minima.
        #This is more obvious to read than mucking about with %
        if tau[-1] < tau[-2] and tau[-1] < tau[0]:
            mins = np.append(mins,-1)
        if tau[0] < tau[-1] and tau[0] < tau[1]:
            mins = np.append(mins,0)
        assert np.size(mins) > 0
        minn = np.max(np.append(mins[np.where(mins < midpt)], 0))
        maxx = np.min(np.append(mins[np.where(mins > midpt)], np.size(tau)))
        assert minn < midpt and maxx > midpt
        diff = np.sum(((tau - gauss))[minn:maxx]**2)
        return diff

    def profile_multiple(self, stddevs, means, amplitudes):
        """Function to generate a profile consisting of multiple peaks, so we can re-fit one last time."""
        gauss = np.sum([self.profile(stddev, mean, amplitude) for (stddev, mean, amplitude) in zip(stddevs, means, amplitudes)],axis=0)
        return gauss

    def fun_min_multiple(self, inputs):
        """Helper function for scipy.optimise that trys to fit the whole (original) spectrum at once
        with the detected peaks."""
        third = int(np.size(inputs)/3)
        stddevs = inputs[0:third]
        means = inputs[third:2*third]
        amplitudes = inputs[2*third:]
        gauss=self.profile_multiple(stddevs, means, amplitudes)
        return np.sum((self.tau - gauss)**2)

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
        mean = mean % np.max(self.wavelengths)
        midpt = int(round(np.size(self.wavelengths)/2))
        T0 = (self.wavelengths - midpt*self.dvbin)/np.abs(stddev)
        aa = self.voigt_fac/np.abs(stddev)
        #This is the Fadeeva function
        fadeeva = np.real(scipy.special.wofz(T0+1j*aa))
        norm = np.real(scipy.special.wofz(1j*aa))
        profile = np.roll(fadeeva/norm*amplitude, int(round(mean/self.dvbin))-midpt)
        assert (np.argmax(profile) - mean/self.dvbin) <= 1
        return profile

    def gaussian_profile(self, stddev,mean,amplitude):
        """A Gaussian profile shape"""
        midpt = int(np.size(self.wavelengths)/2)
        gauss = np.exp(-(self.wavelengths-midpt)**2/stddev**2/2)*amplitude
        profile = np.roll(gauss,int(mean/self.dvbin) - midpt)
        assert (np.argmax(profile) - mean/self.dvbin) <= 1
        return profile

    def get_fitted_profile(self):
        """Helper function to return the fitted profile shape"""
        return (self.wavelengths, self.profile_multiple(self.stddev_new, self.means_new, self.amplitudes_new))

    def get_b_params(self):
        """Helper function to return the doppler b (in km/s)"""
        return self.stddev_new

    def get_positions(self):
        """Helper function to return the wavelength (offset from a 0 at the start of the box) of each absorber."""
        return np.array(self.means_new)

    def get_column_densities(self):
        """Helper function to return the column densities in 1/cm^-2 of each absorber."""
        # amp is a cross-section in cm^2.
        # sqrt(pi)*c / btherm comes from the profile.
        amp = self.sigma_a / math.sqrt(math.pi) * (self.light/self.stddev_new)
        aa = self.voigt_fac/np.abs(self.stddev_new)
        #This is the Fadeeva function
        norm = np.real(scipy.special.wofz(1j*aa))
        #Find amplitude divided by Voigt profile, which is still dimensionless.
        colden = np.array(self.amplitudes_new)/norm/amp
        return colden

def get_voigt_fit_params(taus, dvbin):
    """Helper function to get the Voigt parameters, N_HI and b in a single call."""
    b_params = np.array([])
    n_vals = np.array([])
    for tau_t in taus:
        prof = Profiles(tau_t, dvbin)
        prof.do_fit()
        np.append(n_vals, prof.get_column_densities())
        print("Fit: N=",np.max(prof.get_column_densities()), "b=",prof.get_b_params()[np.argmax(prof.get_column_densities())])
        np.append(b_params, prof.get_b_params())
    return (n_vals, b_params)
