# -*- coding: utf-8 -*-
"""Module for finding the neutral hydrogen in a halo. Each class has a function get_reproc_rhoHI, which returns
the neutral hydrogen density in *physical* atoms / cm^3
    Contains:
        GasProperties: converts reported code outputs to useful physical quantities.
                  get_temp: gets temperature from internal energy
                  get_code_rhoH: gets the neutral density from the code
                  get_reproc_HI - Gets a corrected neutral hydrogen density,
                  so that we are neutral even for star forming gas.
"""

import numpy as np
import scipy.interpolate.interpolate as intp

from . import unitsystem

class GasProperties(object):
    """Class implementing the neutral fraction ala Rahmati 2012.
    Arguments:
        redshift - redshift at which the data was drawn, to compute the self-shielding correction.
        absnap - AbstractSnapshot instance from which to get the particle data.
        hubble - Hubble parameter.
        fbar - Baryon fraction
        units - UnitSystem instance
        sf_neutral - If True (the default) then gas on the star-forming equation of state is assumed to be neutral.
        Should only be true if used with a Springel-Hernquist star formation model in a version of Gadget/Arepo which incorrectly
        sets the neutral fraction in the star forming gas to less than unity.
        """
    def __init__(self, redshift, absnap, hubble = 0.71, fbar=0.17, units=None, sf_neutral=True):
        if units is not None:
            self.units = units
        else:
            self.units = unitsystem.UnitSystem()
        self.absnap = absnap
        self.f_bar = fbar
        self.redshift = redshift
        self.sf_neutral = sf_neutral
        #Interpolate for opacity and gamma_UVB
        #Opacities for the FG09 UVB from Rahmati 2012.
        #IMPORTANT: The values given for z > 5 are calculated by fitting a power law and extrapolating.
        #Gray power law was: -1.12e-19*(zz-3.5)+2.1e-18 fit to z > 2.
        #gamma_UVB was: -8.66e-14*(zz-3.5)+4.84e-13
        #This is clearly wrong, but this model is equally a poor choice at these redshifts anyway.
        gray_opac = [2.59e-18,2.37e-18,2.27e-18, 2.15e-18, 2.02e-18, 1.94e-18, 1.82e-18, 1.71e-18, 1.60e-18]
        gamma_UVB = [3.99e-14, 3.03e-13, 6e-13, 5.53e-13, 4.31e-13, 3.52e-13, 2.678e-13,  1.81e-13, 9.43e-14]
        zz = [0, 1, 2, 3, 4, 5, 6, 7,8]
        self.redshift_coverage = True
        if redshift > zz[-1]:
            self.redshift_coverage = False
            print("Warning: no self-shielding at z=",redshift)
        else:
            gamma_inter = intp.interp1d(zz,gamma_UVB)
            gray_inter = intp.interp1d(zz,gray_opac)
            self.gray_opac = gray_inter(redshift)
            self.gamma_UVB = gamma_inter(redshift)
        #self.hy_mass = 0.76 # Hydrogen massfrac
        self.gamma=5./3
        #Boltzmann constant (cgs)
        self.boltzmann=1.38066e-16
        self.hubble = hubble
        #Physical density threshold for star formation in H atoms / cm^3
        self.PhysDensThresh = self._get_rho_thresh(hubble)

    def _photo_rate(self, nH, temp):
        """Photoionisation rate as  a function of density from Rahmati 2012, eq. 14.
        Calculates Gamma_{Phot}.
        Inputs: hydrogen density, temperature
            n_H
        The coefficients are their best-fit from appendix A."""
        nSSh = self._self_shield_dens(temp)
        photUVBratio= 0.98*(1+(nH/nSSh)**1.64)**-2.28+0.02*(1+nH/nSSh)**-0.84
        return photUVBratio * self.gamma_UVB

    def _self_shield_dens(self,temp):
        """Calculate the critical self-shielding density. Rahmati 202 eq. 13.
        gray_opac and gamma_UVB are parameters of the UVB used.
        gray_opac is in cm^2 (2.49e-18 is HM01 at z=3)
        gamma_UVB in 1/s (1.16e-12 is HM01 at z=3)
        temp is particle temperature in K
        f_bar is the baryon fraction. 0.17 is roughly 0.045/0.265
        Returns density in atoms/cm^3"""
        T4 = temp/1e4
        G12 = self.gamma_UVB/1e-12
        return 6.73e-3 * (self.gray_opac / 2.49e-18)**(-2./3)*(T4)**0.17*(G12)**(2./3)*(self.f_bar/0.17)**(-1./3)

    def _recomb_rate(self, temp):
        """The recombination rate from Rahmati eq A3, also Hui Gnedin 1997.
        Takes temperature in K, returns rate in cm^3 / s"""
        lamb = 315614./temp
        return 1.269e-13*lamb**1.503 / (1+(lamb/0.522)**0.47)**1.923

    def _neutral_fraction(self, nH, temp):
        """The neutral fraction from Rahmati 2012 eq. A8"""
        alpha_A = self._recomb_rate(temp)
        #A6 from Theuns 98
        LambdaT = 1.17e-10*temp**0.5*np.exp(-157809./temp)/(1+np.sqrt(temp/1e5))
        A = alpha_A + LambdaT
        B = 2*alpha_A + self._photo_rate(nH, temp)/nH + LambdaT
        return (B - np.sqrt(B**2-4*A*alpha_A))/(2*A)

    def get_temp(self,part_type, segment):
        """Compute temperature (in K) from internal energy."""
        return self.absnap.get_temp(part_type, segment=segment)

    def get_code_rhoH(self,part_type, segment):
        """Convert density to physical atoms /cm^3: internal gadget density unit is h^2 (1e10 M_sun) / kpc^3"""
        nH = self.absnap.get_data(part_type, "Density", segment=segment)
        conv = np.float32(self.units.UnitDensity_in_cgs*self.hubble**2/(self.units.protonmass)*(1+self.redshift)**3)
        #Convert to physical
        return nH*conv

    def _code_neutral_fraction(self, part_type, segment):
        """Get the neutral fraction from the code"""
        return self.absnap.get_data(part_type, "NeutralHydrogenFraction", segment=segment)

    def get_reproc_HI(self, part_type, segment):
        """Get a neutral hydrogen *fraction* using values given by Arepo
        which are based on Rahmati 2012 if UVB_SELF_SHIELDING is on.
        Above the star formation density use the Rahmati fitting formula directly,
        as Arepo reports values for the eEOS. """
        nH0 = self._code_neutral_fraction(part_type=part_type, segment=segment)
        if not self.sf_neutral:
            return nH0
        #Above star-formation threshold, we want a neutral fraction which includes
        #explicitly the amount of gas in cold clouds.
        #Ideally we should compute this fraction, and then do
        #  tcool = self.get_tcool(nH,bar)[ind]
        #  print np.median(tcool)
        #  cold_frac = self.star.cold_gas_frac(nH[ind], tcool,self.PhysDensThresh/0.76)
        #  print np.mean(cold_frac)
        #  ssnH0 = (1-cold_frac)*self.neutral_fraction(nH[ind], temp[ind])+ cold_frac*self.neutral_fraction(nH[ind], 1e4)
        #But the cooling time reported by the code is not quite what we want here,
        #because it uses the internal energy reported by U, whereas we really want
        #the cooling time using the energy for the hot phase, at a given density.
        #So just assume that at the threshold all gas is in cold clouds.
        #In reality, this should be about 90% of gas in cold clouds, so
        #we will overpredict the neutral fraction by a small amount.
        density = self.absnap.get_data(part_type, "Density", segment=segment)
        conv = np.float32(self.units.UnitDensity_in_cgs*self.hubble**2/(self.units.protonmass)*(1+self.redshift)**3)
        ind = np.where(density > self.PhysDensThresh/0.76/conv)
        if self.redshift_coverage:
            ssnH0 = self._neutral_fraction(density[ind]*conv, 1e4)
            nH0[ind] = ssnH0
        else:
            nH0[ind] = 1.
        return nH0

    def _get_rho_thresh(self, hubble=0.7,t_0_star=2.27,T_SN=5.73e7,T_c = 1000, A_0=573):
        """
        This function calculates the physical density threshold for star formation, using the gadget model.
        See Springel & Hernquist 2003 (astro-ph/0206393) and
        Nagamine, Springel and Hernquist 2004 (astro-ph/0305409).

        Parameters (of the star formation model):
            hubble - hubble parameter in units of 100 km/s/Mpc
            t_0_star - star formation timescale at threshold density
                 - (MaxSfrTimescale) 1.5 in internal time units ( 1 itu ~ 0.97 Gyr/h)
            T_SN - Temperature of the supernova in K- 10^8 K SH03. (TempSupernova) Used to calculate u_SN
            T_c  - Temperature of the cold clouds in K- 10^3 K SH03. (TempClouds) Used to calculate u_c.
            A_0  - Supernova evaporation parameter (FactorEVP = 1000).
            WARNING: cooling time for the cloud is hard-coded, as is rescaled beta.
            This uses values from the default GFM parameters. Do not try to change them!
        Returns:
                rho_thresh in units of H atoms/cm^3
        """
        #Some constants and unit systems
        #Internal velocity unit : 1 km/s in cm/s
        UnitTime_in_s = (self.units.UnitLength_in_cm/self.units.UnitVelocity_in_cm_per_s)
        hy_mass = 0.76 # Primordial Hydrogen massfrac

        #Supernova timescale in s
        t_0_star=t_0_star*UnitTime_in_s/hubble # Now in s

        #This is the value output by an IMF rescaling, if GFM_EVOLUTION is on.
        beta = 0.264089

        #u_c - thermal energy in the cold gas.
        meanweight = 4 / (1 + 3 * hy_mass)          #Assuming neutral gas for u_c
        u_c =  1. / meanweight * (1.0 / (self.gamma-1)) * (self.boltzmann / self.units.protonmass) *T_c

        #SN energy: u_SN = (1-beta)/beta epsilon_SN
        meanweight = 4 / (8 - 5 * (1 - hy_mass))    #Assuming FULL ionization for u_H
        u_SN =  1. / meanweight * (1.0 / (self.gamma -1)) * (self.boltzmann / self.units.protonmass) * T_SN

        #This is the hard-coded "very high density" at which the cooling rate is computed
        #in the Gadget 2-phase SFR model.
        rhoinf = 277.476 * self.units.UnitDensity_in_cgs*self.hubble**2
        #The cooling time for cold clouds, as computed inside the Gadget star-forming model.
        tcool = 4.64419e-10 * UnitTime_in_s/self.hubble

        u_h = u_SN / A_0

        #u_4 - thermal energy at 10^4K
        meanweight = 4 / (8 - 5 * (1 - hy_mass))    #Assuming FULL ionization for u_H
        u_4 =  1. / meanweight * (1.0 / (self.gamma-1)) * (self.boltzmann / self.units.protonmass) *1e4
        coolrate = u_h / tcool / rhoinf

        x = (u_h - u_4) / (u_h - u_c)
        physdens =  x / (1 - x)**2 * (beta * u_SN - (1 -beta) * u_c) /(t_0_star * coolrate)
        return physdens / self.units.protonmass *hy_mass
