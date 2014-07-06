# -*- coding: utf-8 -*-
"""This module stores routines specific to velocity width analysis of fake spectra
"""

import numpy as np
import math
import spectra as ss

class VWSpectra(ss.Spectra):
    """"Extends the spectra class with velocity width functions."""

    def find_absorber_width(self, elem, ion, chunk = 20, minwidth=None):
        """
           Find the region in velocity space considered to be an absorber for each spectrum.
           This is defined to be the maximum of 1000 km/s and the region over which there is "significant"
           absorption in the strongest line for this ion, where strongest is the line with the largest
           cross-section, ie, greatest lambda * fosc.
           elem, ion - ion to look at

           This line will be highly saturated, so consider significant absorption as F < 3/snr,
           or F < 0.15 for no noise (and an assumed SNR of 20).

           Returns the low and high indices of absorption, and the offset for the maximal absorption.
        """
        if minwidth == None:
            minwidth = self.minwidth
        try:
            return self.absorber_width[(elem, ion, minwidth)]
        except KeyError:
            pass
        if self.snr > 0:
            thresh = - np.log(1-4./self.snr)
        else:
            thresh = -np.log(1-0.15)
        lines = self.lines[(elem,ion)]
        strength = [ll.fosc_X*ll.lambda_X for ll in lines.values()]
        ind = np.where(strength == np.max(strength))[0][0]
        #Lines are indexed by wavelength
        strlam = int(lines.values()[ind].lambda_X)
        #Absorption in a strong line: eg, SiII1260.
        strong = self.get_tau(elem, ion, strlam)
        (offset, roll) = self._get_rolled_spectra(strong)
        #Minimum
        if minwidth > 0 and minwidth < self.nbins/2:
            low  = int(self.nbins/2-minwidth/self.dvbin)*np.ones(self.NumLos, dtype=np.int)
            high = int(self.nbins/2+minwidth/self.dvbin)*np.ones(self.NumLos, dtype=np.int)
        else:
            low = np.zeros(self.NumLos, dtype=np.int)
            high = self.nbins*np.ones(self.NumLos, dtype=np.int)
        for ii in xrange(self.NumLos):
            #First expand the search area in case there is absorption at the edges.
            for i in xrange(low[ii],0,-chunk):
                if not np.any(roll[ii,i:(i+chunk)] > thresh):
                    low[ii] = i
                    break
            #Where is there no absorption rightwards of the peak?
            for i in xrange(high[ii],self.nbins,chunk):
                if not np.any(roll[ii,i:(i+chunk)] > thresh):
                    high[ii] = i+chunk
                    break
            #Shrink to width which has some absorption
            ind = np.where(roll[ii][low[ii]:high[ii]] > thresh)[0]
            if np.size(ind) != 0:
                oldlow = low[ii]
                low[ii] = np.max((ind[0]+oldlow,0))
                high[ii] = np.min((ind[-1]+oldlow+chunk,self.nbins))
        self.absorber_width[(elem, ion, minwidth)] = (low, high, offset)
        return (low, high, offset)

    def _eq_width_from_colden(self, col_den, elem = "H", ion = 1, line = 1215):
        """Find the equivalent width of the line from the column density,
           assuming we are in the damping wing regime. Default line is Lyman-alpha.
           Returns width in km/s.
        """
        #line properties
        line = self.lines[(elem,ion)][line]
        #Convert wavelength to cm
        lambdacgs = line.lambda_X*1e-8
        #Tompson cross-section for the electron
        sigma_t = 6.652458558e-25
        #Line cross-section
        sigma_a = np.sqrt(3*math.pi*sigma_t/8.)*lambdacgs*line.fosc_X
        #In the damping wing, W ~ sqrt(N).
        #Coefficients come from setting tau = 1, and expanding the Voigt function used
        #in Tepper-Garcia 2006 where exp(-x^2) ~ 0 (ie, far from the line center)
        #then simplifying a bit
        width = np.sqrt(line.gamma_X*lambdacgs*self.light*col_den*sigma_a)/math.pi
        #Convert from cm/s to km/s
        return width/1e5

    def get_observer_tau(self, elem, ion, number=-1, force_recompute=False, noise=True):
        """Get the optical depth for a particular element out of:
           (He, C, N, O, Ne, Mg, Si, Fe)
           and some ion number, choosing the line which causes the maximum optical depth to be closest to unity.
        """
        try:
            if force_recompute:
                raise KeyError
            self._really_load_array((elem, ion), self.tau_obs, "tau_obs")
            ntau = self.tau_obs[(elem, ion)]
        except KeyError:
            #Compute tau for each line
            nlines = len(self.lines[(elem,ion)])
            tau = np.empty([nlines, self.NumLos,self.nbins])
            pos = {}
            vel = {}
            elem_den = {}
            temp = {}
            hh = {}
            amumass = {}
            for ff in self.files:
                (pos[ff], vel[ff], elem_den[ff], temp[ff], hh[ff], amumass[ff]) = self._read_particle_data(ff, elem, ion,True)

            for ll in xrange(nlines):
                line = (self.lines[(elem,ion)].values())[ll]
                for ff in self.files:
                    if amumass[ff] != False:
                        tau_loc = self._do_interpolation_work(pos[ff], vel[ff], elem_den[ff], temp[ff], hh[ff], amumass[ff], line, True)
                        tau[ll,:,:] += tau_loc
                        del tau_loc
            #Maximum tau in each spectra with each line,
            #after convolving with a Gaussian for instrumental broadening.
            maxtaus = np.max(self.res_corr(tau,self.spec_res), axis=-1)
            #Array for line indices
            ntau = np.empty([self.NumLos, self.nbins])
            #Use the maximum unsaturated optical depth
            for ii in xrange(self.NumLos):
                # we want unsaturated lines, defined as those with tau < 3
                #which is the maximum tau in the sample of Neeleman 2013
                #Also use lines with some absorption: tau > 0.1, roughly twice noise level.
                ind = np.where(np.logical_and(maxtaus[:,ii] < 3, maxtaus[:,ii] > 0.1))
                if np.size(ind) > 0:
                    line = np.where(maxtaus[:,ii] == np.max(maxtaus[ind,ii]))
                else:
                    #We have no lines in the desired region: here use something slightly saturated.
                    #In reality the observers will use a different ion
                    ind2 = np.where(maxtaus[:,ii] > 0.1)
                    if np.size(ind2) > 0:
                        line = np.where(maxtaus[:,ii] == np.min(maxtaus[ind2,ii]))
                    else:
                        #We have no observable lines: this spectra are metal-poor
                        #and will be filtered anyway.
                        line = np.where(maxtaus[:,ii] == np.max(maxtaus[:,ii]))
                if np.size(line) > 1:
                    line = (line[0][0],)
                ntau[ii,:] = tau[line,ii,:]
            self.tau_obs[(elem, ion)] = ntau
        if number >= 0:
            ntau = ntau[number,:]
        # Convolve lines by a Gaussian filter of the resolution of the spectrograph.
        ntau = self.res_corr(ntau, self.spec_res)
        #Add noise
        if noise and self.snr > 0:
            ntau = self.add_noise(self.snr, ntau, number)
        return ntau

    def vel_width(self, elem, ion):
        """
           Find the velocity width of an ion.
           This is the width in velocity space containing 90% of the optical depth
           over the absorber.

           elem - element to look at
           ion - ionisation state of this element.
        """
        try:
            return self.vel_widths[(elem, ion)]
        except KeyError:
            tau = self.get_observer_tau(elem, ion)
            (low, high, offset) = self.find_absorber_width(elem, ion)
            #  Size of a single velocity bin
            vel_width = np.zeros(np.shape(tau)[0])
            #deal with periodicity by making sure the deepest point is in the middle
            for ll in np.arange(0, np.shape(tau)[0]):
                tau_l = np.roll(tau[ll,:],offset[ll])[low[ll]:high[ll]]
                (nnlow, nnhigh) = self._vel_width_bound(tau_l)
                vel_width[ll] = self.dvbin*(nnhigh-nnlow)
            #Return the width
            self.vel_widths[(elem, ion)] = vel_width
            return self.vel_widths[(elem, ion)]

    def _vel_width_bound(self, tau):
        """Find the 0.05 and 0.95 bounds of the integrated optical depth"""
        #Zero everything less than 1 sigma significant
        cum_tau = np.cumsum(tau)
        #Use spline interpolation to find the edge of the bins.
        tdiff = cum_tau - 0.95*cum_tau[-1]
        high = np.where(tdiff >= 0)[0][0]
        tdiff = cum_tau - 0.05*cum_tau[-1]
        low = np.where(tdiff >= 0)[0][0]
        return (low, high)

    def _vel_median(self, tau):
        """Find the median point of the integrated optical depth"""
        cum_tau = np.cumsum(tau)
        #Use spline interpolation to find the edge of the bins.
        tdiff = cum_tau - 0.5*cum_tau[-1]
        high = np.where(tdiff >= 0)[0][0]
        return high

    def vel_mean_median(self, elem, ion):
        """Find the difference between the mean velocity and the median velocity.
           The mean velocity is the point halfway across the extent of the velocity width.
           The median velocity is v(tau = tot_tau /2)
           """
        tau = self.get_observer_tau(elem, ion)
        (low, high, offset) = self.find_absorber_width(elem, ion)
        mean_median = np.zeros(np.shape(tau)[0])
        #Deal with periodicity by making sure the deepest point is in the middle
        for ll in xrange(np.shape(tau)[0]):
            tau_l = np.roll(tau[ll,:],offset[ll])[low[ll]:high[ll]]
            (nnlow, nnhigh) = self._vel_width_bound(tau_l)
            vel_median = self._vel_median(tau_l)
            vmean = (nnlow+nnhigh)/2.
            mean_median[ll] = np.abs(vmean - vel_median)/((nnhigh-nnlow)/2.)
        #Return the width
        return mean_median

    def vel_peak(self, elem, ion):
        """
           Find the f_peak statistic for spectra in an ion.
           f_peak = (vel_peak - vel_mean) / (v_90/2)
        """
        tau = self.get_observer_tau(elem, ion)
        (low, high, offset) = self.find_absorber_width(elem, ion)
        peak = np.zeros(np.shape(tau)[0])
        #Deal with periodicity by making sure the deepest point is in the middle
        for ll in xrange(np.shape(tau)[0]):
            tau_l = np.roll(tau[ll,:],offset[ll])[low[ll]:high[ll]]
            peak[ll] = self._vel_peak_tau(tau_l)
        #Return the width
        return peak

    def _vel_peak_tau(self,tau_l):
        """Helper function for a single spectrum to compute v_peak"""
        (nnlow, nnhigh) = self._vel_width_bound(tau_l)
        vmax = np.where(tau_l == np.max(tau_l))[0][0]
        vmean = (nnlow+nnhigh)/2.
        peak = np.abs(vmax - vmean)/((nnhigh-nnlow)/2.)
        return peak

    def vel_width_hist(self, elem, ion, dv=0.1):
        """
        Compute a histogram of the velocity widths of our spectra, with the purpose of
        comparing to the data of Prochaska 2008.

        Note this does not match Pontzen 2008, who multiply by the DLA fraction (0.065) obtained from the cddf.

        So we have f(N) = d n/ dv
        and n(N) = number of absorbers per sightline in this velocity bin.
        Note f(N) has dimensions of s/km, because v has units of km/s.

        Parameters:
            elem - element to use
            line - line to use (the components of this line must be pre-computed and stored in self.metals)
            dv - bin spacing
            met_cut - Discard spectra whose maximal metal column density is below this level.
                      Removes unobservable systems.

        Returns:
            (v, f_table) - v (binned in log) and corresponding f(N)
        """
        return self._vel_stat_hist(elem, ion, dv, self.vel_width, log=True)

    def f_meanmedian_hist(self, elem, ion, dv=0.1):
        """
        Compute a histogram of the mean median statistic of our spectra, the difference in
        units of the velocity width between the mean velocity and median velocity of
        the absorber.

        For arguments see vel_width_hist.
        """
        return self._vel_stat_hist(elem, ion, dv, self.vel_mean_median, log=False)

    def f_peak_hist(self, elem, ion, dv=0.1):
        """
        Compute a histogram of the peak statistic of our spectra, the difference in
        units of the velocity width between the largest peak velocity and the mean velocity of
        the absorber.

        For arguments see vel_width_hist.
        """
        return self._vel_stat_hist(elem, ion, dv, self.vel_peak, log=False)

    def get_separated(self, elem="Si", ion = 2, thresh = 1e-1, mindist=0):
        """
        Find spectra with more than a single density peak.
        Threshold is as a percentage of the maximum value.
        mindist is in km/s
        """
        dist = int(mindist/self.dvbin)
        ind = self.get_filt(elem, ion)
        rho = self.get_col_density(elem, ion)[ind]
        H1_den = self.get_col_density("H", 1)[ind]
        seps = np.zeros(np.size(ind[0]), dtype=np.bool)
        lls = 0
        dla = 0
        none = 0
        #deal with periodicity by making sure the deepest point is in the middle
        for ll in xrange(np.size(ind[0])):
            rho_l = rho[ll,:]
            H1_l = H1_den[ll,:]
            lsep = ss.combine_regions(rho_l > thresh*np.max(rho_l), dist)
            seps[ll] = (np.shape(lsep)[0] > 1)
            if seps[ll] == False:
                continue
            m_H1_colden = np.array([np.sum(H1_l[lsep[jj,0]:lsep[jj,1]]) for jj in xrange(np.shape(lsep)[0])])
            #All DLAs
            if np.all(m_H1_colden > 10**(20.3)):
                dla += 1
            #Some LLS
            elif np.all(m_H1_colden > 10**(17.)):
                lls += 1
            else:
                none += 1
        tot = dla + lls + none
        print "Fraction DLA: ",1.*dla/tot," Fraction LLS: ",1.*lls/tot," fraction less: ",1.*none/tot
        return seps

