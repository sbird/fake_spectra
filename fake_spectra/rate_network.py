# -*- coding: utf-8 -*-
"""A rate network for neutral hydrogen following
Katz, Weinberg & Hernquist 1996, eq. 28-32."""
import os.path
import math
import numpy as np
import scipy.interpolate as interp

class RateNetwork(object):
    """A rate network for neutral hydrogen following
    Katz, Weinberg & Hernquist 1996, astro-ph/9509107, eq. 28-32.

    Most internal methods are CamelCapitalized and follow a convention that
    they are named like the process and then the ion they refer to.
    eg:
        CollisionalExciteHe0 is the neutral Helium collisional excitation rate.
        RecombHp is the recombination rate for ionized hydrogen.

    Externally useful methods (the API) are named like get_*.
    These are:
        get_temp() - gets the temperature from the density and internal energy.
        get_cooling_rate() - gets the total cooling rate from density and internal energy.
        get_neutral_fraction() - gets the neutral fraction from the rate network given density and internal energy.
    Two useful helper functions:
        get_equilib_ne() - gets the equilibrium electron density.
        get_ne_by_nh() - gets the above, divided by the hydrogen density (Gadget reports this as ElectronAbundance).

    Constructor arguments:
        redshift - the redshift at which to evaluate the cooling. Affects the photoionization rate,
                   the Inverse Compton cooling and the self shielding threshold.
        photo_factor - Factor by which to multiply the UVB amplitude.
        f_bar - Baryon fraction. Omega_b / Omega_cdm.
        converge - Tolerance to which the rate network should be converged.
        selfshield - Flag to enable self-shielding following Rahmati 2013
        cool - which cooling rate coefficient table to use.
               Supported are: KWH (original Gadget rates)
                              Nyx (rates used in Nyx (Lukic 2015))
                              Sherwood (rates used in Sherwood simulations (Bolton 2017))
              Default is Sherwood
        recomb - which recombination rate table to use.
                 Supported are: C92 (Cen 1992, the Gadget default)
                                V96 (Verner & Ferland 1996, more accurate rates).
                                B06 (Badnell 2006 rates, current cloudy defaults. Very similar to V96).
        collisional - Flag to enable collisional ionizations.
        treecool_file - File to read a UV background from. Matches format used by Gadget.
    """
    def __init__(self,redshift, photo_factor = 1., f_bar = 0.17, converge = 1e-7, selfshield=True, cool="Sherwood", recomb="V96", collisional=True, treecool_file="data/TREECOOL_ep_2018p"):
        if recomb == "V96":
            self.recomb = RecombRatesVerner96()
        elif recomb == "B06":
            self.recomb = RecombRatesBadnell()
        else:
            self.recomb = RecombRatesCen92()
        self.photo = PhotoRates(treecool_file=treecool_file)
        self.photo_factor = photo_factor
        self.f_bar = f_bar
        if cool == "KWH":
            self.cool = CoolingRatesKWH92()
        elif cool == "Sherwood":
            self.cool = CoolingRatesSherwood()
        elif cool == "Nyx":
            self.cool = CoolingRatesNyx()
        else:
            raise ValueError("Not supported")
        #Extra helium reionization photoheating model
        self.hub = 0.7
        self.he_thresh = 10
        self.he_amp = 1
        self.he_exp = 0
        self.he_model_on = False
        #proton mass in g
        self.protonmass = 1.67262178e-24
        self.redshift = redshift
        self.converge = converge
        self.selfshield = selfshield
        self.collisional = collisional
        zz = [0, 1, 2, 3, 4, 5]
        #Tables for the self-shielding correction.
        #These are not well-measured for z > 5, so we keep the ss density constant.
        gray_opac = [2.59e-18,2.37e-18,2.27e-18, 2.15e-18, 2.02e-18, 1.94e-18]
        self.Gray_ss = interp.interp1d(zz, gray_opac, bounds_error=False, fill_value = (gray_opac[0], gray_opac[-1]))

    def get_temp(self, density, ienergy, helium=0.24):
        """Get the equilibrium temperature at given internal energy.
           Density is gas density in (physical) protons/cm^3
           Internal energy is in (km/s)^2 (internal gadget units) == 10^10 ergs/g.
        helium is a mass fraction"""
        ne = self.get_equilib_ne(density, ienergy, helium)
        nh = density * (1-helium)
        return self._get_temp(ne/nh, ienergy, helium)

    def get_cooling_rate(self, density, ienergy, helium=0.24, photoheating=False):
        """Get the total cooling rate for a temperature and density. Negative means heating.
           Density is gas density in (physical) protons/cm^3
           Internal energy is in (km/s)^2 (internal gadget units) == 10^10 ergs/g.
           Returns net heating/cooling rate in erg/s/g.
        """
        ne = self.get_equilib_ne(density, ienergy, helium)
        nh = density * (1-helium)
        nebynh = ne / nh
        temp = self._get_temp(nebynh, ienergy, helium)
        nH0 = self._nH0(nh, temp, nebynh)
        nHp = self._nHp(nh, temp, nebynh)
        yy = helium / 4 / (1 - helium)
        nHe0 = self._nHe0(nh, temp, nebynh) * yy
        nHep = self._nHep(nh, temp, nebynh) * yy
        nHepp = self._nHepp(nh, temp, nebynh) * yy
        #This is the collisional excitation and ionisation rate.
        LambdaCollis = nebynh * (self.cool.CollisionalH0(temp) * nH0 +
                                self.cool.CollisionalHe0(temp) * nHe0 +
                                self.cool.CollisionalHeP(temp) * nHep)
        LambdaRecomb = nebynh * (self.cool.RecombHp(temp) * nHp +
                                self.cool.RecombHeP(temp) * nHep +
                                self.cool.RecombHePP(temp) * nHepp)
        LambdaFF = nebynh * (self.cool.FreeFree(temp, 1)*(nHp + nHep) + self.cool.FreeFree(temp, 2)*nHepp)
        LambdaCmptn = nebynh**2 * self.cool.InverseCompton(temp, self.redshift)
        Lambda = LambdaCollis + LambdaRecomb + LambdaFF + LambdaCmptn
        Heating = 0
        if photoheating:
            Heating = nH0 * self.photo.epsH0(self.redshift)/nh
            Heating += nHe0 * self.photo.epsHe0(self.redshift)/nh
            Heating += nHep * self.photo.epsHep(self.redshift)/nh
            Heating *= self.photo_factor
            if self.he_model_on:
                Heating *= self._he_reion_factor(density)
        #print("Heat = %g Lambda = %g LC = %g LR = %g LFF = %g LCmptn = %g, ne = %g, nH0 = %g, nHp = %g, nHe0 = %g, nHep = %g, nHepp = %g, nh=%g, temp=%g, ienergy=%g" % (Heating, Lambda, LambdaCollis, LambdaRecomb, LambdaFF, LambdaCmptn, ne/nh, nH0, nHp, nHe0, nHep, nHepp, nh, temp, ienergy))
        LambdaNet = Lambda - Heating
        # LambdaNet in erg/s cm^3, Density in protons/cm^3, PROTONMASS in protons/g.
        # Convert to erg/s/g*/
        return LambdaNet * (1 - helium)**2 * density / self.protonmass

    def get_equilib_ne(self, density, ienergy,helium=0.24):
        """Solve the system of equations for photo-ionisation equilibrium,
        starting with ne = nH and continuing until convergence.
        Density is gas density in (physical) protons/cm^3
        Internal energy is in (km/s)^2 (internal gadget units) == 10^10 ergs/g.
        helium is a mass fraction.
        """
        #Get hydrogen number density
        nh = density * (1-helium)
        def rooted(nebynh, nh, ienergy):
            return self._nebynh(nh, self._get_temp(nebynh, ienergy, helium=helium), nebynh, helium=helium)
        nebynh = fixed_point(rooted, nh, args=(nh, ienergy),xtol=self.converge)
        assert np.all(np.abs(rooted(nebynh, nh, ienergy) - nebynh) < self.converge)
        return nebynh * nh

    def get_ne_by_nh(self, density, ienergy, helium=0.24):
        """Same as above, but get electrons per proton."""
        return self.get_equilib_ne(density, ienergy, helium)/(density*(1-helium))

    def get_neutral_fraction(self, density, ienergy, helium=0.24):
        """Get the neutral hydrogen fraction at a given temperature and density.
        Density is gas density in (physical) protons/cm^3
        Internal energy is in (km/s)^2 (internal gadget units) == 10^10 ergs/g.
        helium is a mass fraction.
        """
        ne = self.get_equilib_ne(density, ienergy, helium=helium)
        nh = density * (1-helium)
        temp = self._get_temp(ne/nh, ienergy, helium)
        return self._nH0(nh, temp, ne/nh)


    def _photofac(self, nebynh, nh, temp):
        """Get the photoionization correction factor divided by the electron density.
        Adjusted to work when ne ~ 0.
        Arguments:
            ne - electron abundance per atom
            nh - proton abundance
            temp - temperature"""
        if np.size(nebynh) > 1:
            photofac = np.zeros_like(nebynh)
            ii = np.where(nebynh > 1e-40)
            photofac[ii] = self.photo_factor*self._self_shield_corr(nh[ii], temp[ii])/(nebynh[ii]*nh[ii])
        else:
            photofac = 0
            if nebynh > 1e-40:
                photofac = self.photo_factor*self._self_shield_corr(nh, temp)/(nebynh*nh)
        return photofac

    def _nH0(self, nh, temp, nebynh):
        """The neutral hydrogen number density per hydrogen atom. Eq. 33 of KWH."""
        alphaHp = self.recomb.alphaHp(temp)
        GammaeH0 = self.collisional * self.recomb.GammaeH0(temp)
        photorate = self.photo.gH0(self.redshift)* self._photofac(nebynh, nh, temp)
        return alphaHp/ (alphaHp + GammaeH0 + photorate)

    def _nHp(self, nh, temp, nebynh):
        """The ionised hydrogen number density. Eq. 34 of KWH."""
        return 1 - self._nH0(nh, temp, nebynh)

    def _nHep(self, nh, temp, nebynh):
        """The ionised helium number density, divided by the helium number fraction. Eq. 35 of KWH."""
        alphaHep = self.recomb.alphaHep(temp) + self.recomb.alphad(temp)
        alphaHepp = self.recomb.alphaHepp(temp)
        photofac = self._photofac(nebynh, nh, temp)
        GammaHe0 = self.collisional * self.recomb.GammaeHe0(temp) + self.photo.gHe0(self.redshift)*photofac
        GammaHep = self.collisional * self.recomb.GammaeHep(temp) + self.photo.gHep(self.redshift)*photofac
        return 1. / (1 + alphaHep / GammaHe0 + GammaHep/alphaHepp)

    def _nHe0(self, nh, temp, nebynh):
        """The neutral helium number density, divided by the helium number fraction. Eq. 36 of KWH."""
        alphaHep = self.recomb.alphaHep(temp) + self.recomb.alphad(temp)
        photofac = self._photofac(nebynh, nh, temp)
        GammaHe0 = self.collisional * self.recomb.GammaeHe0(temp) + self.photo.gHe0(self.redshift)*photofac
        return self._nHep(nh, temp, nebynh) * alphaHep / GammaHe0

    def _nHepp(self, nh, temp, nebynh):
        """The doubly ionised helium number density, divided by the helium number fraction. Eq. 37 of KWH."""
        photofac = self._photofac(nebynh, nh, temp)
        GammaHep = self.collisional * self.recomb.GammaeHep(temp) + self.photo.gHep(self.redshift)*photofac
        alphaHepp = self.recomb.alphaHepp(temp)
        return self._nHep(nh, temp, nebynh) * GammaHep / alphaHepp

    def _nebynh(self, nh, temp, nebynh, helium=0.24):
        """The electron number density per hydrogen atom. Eq. 38 of KWH."""
        yy = helium / 4 / (1 - helium)
        return self._nHp(nh, temp, nebynh) + yy * self._nHep(nh, temp, nebynh) + 2* yy * self._nHepp(nh, temp, nebynh)

    def _self_shield_corr(self, nh, temp):
        """Photoionisation rate as  a function of density from Rahmati 2012, eq. 14.
        Calculates Gamma_{Phot} / Gamma_{UVB}.
        Inputs: hydrogen density, temperature
            n_H
        The coefficients are their best-fit from appendix A."""
        if not self.selfshield:
            return np.ones_like(nh)
        nSSh = 1.003*self._self_shield_dens(self.redshift, temp)
        return 0.98*(1+(nh/nSSh)**1.64)**-2.28+0.02*(1+nh/nSSh)**-0.84

    def _self_shield_dens(self,redshift, temp):
        """Calculate the critical self-shielding density. Rahmati 202 eq. 13.
        gray_opac is a parameter of the UVB used.
        gray_opac is in cm^2 (2.49e-18 is HM01 at z=3)
        temp is particle temperature in K
        f_bar is the baryon fraction. 0.17 is roughly 0.045/0.265
        Returns density in atoms/cm^3"""
        T4 = temp/1e4
        G12 = self.photo.gH0(redshift)/1e-12
        return 6.73e-3 * (self.Gray_ss(redshift) / 2.49e-18)**(-2./3)*(T4)**0.17*(G12)**(2./3)*(self.f_bar/0.17)**(-1./3)

    def _he_reion_factor(self, density):
        """Compute a density dependent correction factor to the heating rate which can model the effect of helium reionization.
        Argument: Gas density in protons/cm^3."""
        #Newton's constant (cgs units)
        gravity = 6.672e-8
        #100 km/s/Mpc in h/sec
        hubble = 3.2407789e-18
        omegab = 0.0483
        atime = 1/(1+self.redshift)
        rhoc = 3 * (self.hub* hubble)**2 /(8* math.pi * gravity)
        overden = self.protonmass * density /(omegab * rhoc * atime**(-3))
        if overden >= self.he_thresh:
            overden = self.he_thresh
        return self.he_amp * overden**self.he_exp

    def _get_temp(self, nebynh, ienergy, helium=0.24):
        """Compute temperature (in K) from internal energy and electron density.
           Uses: internal energy in (km/s)^2, internal gadget units.
                 electron abundance per H atom (ne/nH)
                 hydrogen mass fraction (0.76)
           Internal energy is in (km/s)^2 (internal gadget units) == 10^10 ergs/g.
           The code takes internal Gadget units
           Factor to convert U (J/kg) to T (K) : U = N k T / (γ - 1)
           T = U (γ-1) μ m_P / k_B
           where k_B is the Boltzmann constant
           γ is 5/3, the perfect gas constant
           m_P is the proton mass

           μ = 1 / (mean no. molecules per unit atomic weight)
             = 1 / (X + Y /4 + E)
             where E = Ne * X, and Y = (1-X).
             Can neglect metals as they are heavy.
             Leading contribution is from electrons, which is already included
             [+ Z / (12->16)] from metal species
             [+ Z/16*4 ] for OIV from electrons."""
        #convert U (J/kg) to T (K) : U = N k T / (γ - 1)
        #T = U (γ-1) μ m_P / k_B
        #where k_B is the Boltzmann constant
        #γ is 5/3, the perfect gas constant
        #m_P is the proton mass
        #μ is 1 / (mean no. molecules per unit atomic weight) calculated in loop.
        #Internal energy units are 10^-10 erg/g
        hy_mass = 1 - helium
        meanweight = 4./(1 + (3 + 4 * nebynh) * hy_mass)
        muienergy = meanweight *ienergy*1e10
        #Boltzmann constant (cgs)
        boltzmann=1.38066e-16
        gamma=5./3
        #So for T in K, boltzmann in erg/K, internal energy has units of erg/g
        temp = (gamma-1) * self.protonmass / boltzmann * muienergy
        assert np.all(temp >= 0)
        return temp

class RecombRatesCen92(object):
    """Recombination rates and collisional ionization rates, as a function of temperature.
    This is taken from KWH 06, astro-ph/9509107, Table 2, based on Cen 1992.
    Illustris uses these rates."""
    def alphaHp(self,temp):
        """Recombination rate for H+, ionized hydrogen, in cm^3/s.
        Temp in K."""
        return 8.4e-11 / np.sqrt(temp) / np.power(temp/1000, 0.2) / (1+ np.power(temp/1e6, 0.7))

    def alphaHep(self,temp):
        """Recombination rate for He+, ionized helium, in cm^3/s.
        Temp in K."""
        return 1.5e-10 / np.power(temp,0.6353)

    def alphad(self, temp):
        """Recombination rate for dielectronic recombination, in cm^3/s.
        Temp in K."""
        return 1.9e-3 / np.power(temp,1.5) * np.exp(-4.7e5/temp)*(1+0.3*np.exp(-9.4e4/temp))

    def alphaHepp(self, temp):
        """Recombination rate for doubly ionized helium, in cm^3/s.
        Temp in K."""
        return 4 * self.alphaHp(temp)

    def GammaeH0(self,temp):
        """Collisional ionization rate for H0 in cm^3/s. Temp in K"""
        return 5.85e-11 * np.sqrt(temp) * np.exp(-157809.1/temp) / (1+ np.sqrt(temp/1e5))

    def GammaeHe0(self,temp):
        """Collisional ionization rate for H0 in cm^3/s. Temp in K"""
        return 2.38e-11 * np.sqrt(temp) * np.exp(-285335.4/temp) / (1+ np.sqrt(temp/1e5))

    def GammaeHep(self,temp):
        """Collisional ionization rate for H0 in cm^3/s. Temp in K"""
        return 5.68e-12 * np.sqrt(temp) * np.exp(-631515.0/temp) / (1+ np.sqrt(temp/1e5))

class RecombRatesVerner96(object):
    """Recombination rates and collisional ionization rates, as a function of temperature.
     Recombination rates are the fit from Verner & Ferland 1996 (astro-ph/9509083).
     Collisional rates are the fit from Voronov 1997 (http://www.sciencedirect.com/science/article/pii/S0092640X97907324).

     In a very photoionised medium this changes the neutral hydrogen abundance by approximately 10% compared to Cen 1992.
     These rates are those used by Nyx.
    """
    def _Verner96Fit(self, temp, aa, bb, temp0, temp1):
        """Formula used as a fitting function in Verner & Ferland 1996 (astro-ph/9509083)."""
        sqrttt0 = np.sqrt(temp/temp0)
        sqrttt1 = np.sqrt(temp/temp1)
        return aa / ( sqrttt0 * (1 + sqrttt0)**(1-bb)*(1+sqrttt1)**(1+bb) )

    def alphaHp(self,temp):
        """Recombination rate for H+, ionized hydrogen, in cm^3/s.
        The V&F 96 fitting formula is accurate to < 1% in the worst case.
        Temp in K."""
        #See line 1 of V&F96 table 1.
        return self._Verner96Fit(temp, aa=7.982e-11, bb=0.748, temp0=3.148, temp1=7.036e+05)

    def alphaHep(self,temp):
        """Recombination rate for He+, ionized helium, in cm^3/s.
        Accurate to ~2% for T < 10^6 and 5% for T< 10^10.
        Temp in K."""
        #VF96 give two rates. The first is more accurate for T < 10^6, the second is valid up to T = 10^10.
        #We use the most accurate allowed. See lines 2 and 3 of Table 1 of VF96.
        lowTfit = self._Verner96Fit(temp, aa=3.294e-11, bb=0.6910, temp0=1.554e+01, temp1=3.676e+07)
        highTfit = self._Verner96Fit(temp, aa=9.356e-10, bb=0.7892, temp0=4.266e-02, temp1=4.677e+06)
        #Note that at 10^6K the two fits differ by ~10%. This may lead one to disbelieve the quoted accuracies!
        #We thus switch over at a slightly lower temperature.
        #The two fits cross at T ~ 3e5K.
        swtmp = 7e5
        deltat = 1e5
        upper = swtmp + deltat
        lower = swtmp - deltat
        #In order to avoid a sharp feature at 10^6 K, we linearly interpolate between the two fits around 10^6 K.
        interpfit = (lowTfit * (upper - temp) + highTfit * (temp - lower))/(2*deltat)
        return (temp < lower)*lowTfit + (temp > upper)*highTfit + (upper > temp)*(temp > lower)*interpfit

    def alphad(self, temp):
        """Recombination rate for dielectronic recombination, in cm^3/s.
        This is the value from Aldrovandi & Pequignot 73, as used in Nyx, Sherwood and Cen 1992.
        It is corrected from the value in Aldrovandi & Pequignot 1973 by Burgess & Tworkowski 1976 (fig1)
        by a factor of 0.65. The exponent is also made slightly more accurate.
        Temp in K."""
        return 1.23e-3 / np.power(temp,1.5) * np.exp(-4.72e5/temp)*(1+0.3*np.exp(-9.4e4/temp))

    def alphaHepp(self, temp):
        """Recombination rate for doubly ionized helium, in cm^3/s. Accurate to 2%.
        Temp in K."""
        #See line 4 of V&F96 table 1.
        return self._Verner96Fit(temp, aa=1.891e-10, bb=0.7524, temp0=9.370, temp1=2.774e6)

    def _Voronov96Fit(self, temp, dE, PP, AA, XX, KK):
        """Fitting function for collisional rates. Eq. 1 of Voronov 1997. Accurate to 10%,
        but data is only accurate to 50%."""
        bolevk = 8.61734e-5 # Boltzmann constant in units of eV/K
        UU = dE / (bolevk * temp)
        return AA * (1 + PP * np.sqrt(UU))/(XX+UU) * UU**KK * np.exp(-UU)

    def GammaeH0(self,temp):
        """Collisional ionization rate for H0 in cm^3/s. Temp in K. Voronov 97, Table 1."""
        return self._Voronov96Fit(temp, 13.6, 0, 0.291e-07, 0.232, 0.39)

    def GammaeHe0(self,temp):
        """Collisional ionization rate for He0 in cm^3/s. Temp in K. Voronov 97, Table 1."""
        return self._Voronov96Fit(temp, 24.6, 0, 0.175e-07, 0.180, 0.35)

    def GammaeHep(self,temp):
        """Collisional ionization rate for HeI in cm^3/s. Temp in K. Voronov 97, Table 1."""
        return self._Voronov96Fit(temp, 54.4, 1, 0.205e-08, 0.265, 0.25)

class RecombRatesBadnell(RecombRatesVerner96):
    """Recombination rates and collisional ionization rates, as a function of temperature.
     Recombination rates are the fit from Badnell's website: http://amdpp.phys.strath.ac.uk/tamoc/RR/#partial.
    """

    def _RecombRateFit_lowcharge_ion(self, temp, aa, bb, cc, temp0, temp1, temp2):
        """Formula used as a fitting function in Verner & Ferland 1996 (astro-ph/9509083)/ See http://amdpp.phys.strath.ac.uk/tamoc/RR/#partial."""
        sqrttt0 = np.sqrt(temp/temp0)
        sqrttt1 = np.sqrt(temp/temp1)
        BB = bb + cc*np.exp(-temp2/temp)
        return aa / ( sqrttt0 * (1 + sqrttt0)**(1-BB)*(1+sqrttt1)**(1+BB) )


    def alphaHp(self,temp):
        """Recombination rate for H+, ionized hydrogen, in cm^3/s.
        Temp in K."""
        #See line 1 of V&F96 table 1.
        return self._Verner96Fit(temp, aa=8.318e-11, bb=0.7472, temp0=2.965, temp1=7.001e5)

    def alphaHep(self,temp):
        """Recombination rate for H+, ionized hydrogen, in cm^3/s.
        Temp in K."""
        #See line 1 of V&F96 table 1.
        return self._Verner96Fit(temp, aa=1.818E-10, bb=0.7492, temp0=10.17, temp1=2.786e6)

    def alphaHepp(self, temp):
        """Recombination rate for doubly ionized helium, in cm^3/s.
        Temp in K."""
        #See line 4 of V&F96 table 1.
        return self._RecombRateFit_lowcharge_ion(temp, aa=5.235E-11, bb=0.6988, cc=0.0829, temp0=7.301, temp1=4.475e6, temp2 = 1.682e5)

class PhotoRates(object):
    """The photoionization rates for a given species.
    Eq. 29 of KWH 96. This is loaded from a TREECOOL table."""
    def __init__(self, treecool_file="data/TREECOOL_ep_2018p"):
        #Format of the treecool table:
        # log_10(1+z), Gamma_HI, Gamma_HeI, Gamma_HeII,  Qdot_HI, Qdot_HeI, Qdot_HeII,
        # where 'Gamma' is the photoionization rate and 'Qdot' is the photoheating rate.
        # The Gamma's are in units of s^-1, and the Qdot's are in units of erg s^-1.
        try:
            data = np.loadtxt(treecool_file)
        except OSError:
            treefile = os.path.join(os.path.dirname(os.path.realpath(__file__)), treecool_file)
            data = np.loadtxt(treefile)
        redshifts = data[:,0]
        photo_rates = data[:,1:4]
        photo_heat = data[:,4:7]
        assert np.shape(redshifts)[0] == np.shape(photo_rates)[0]
        self.Gamma_HI = interp.interp1d(redshifts, np.log10(photo_rates[:,0]))
        self.Gamma_HeI = interp.interp1d(redshifts, np.log10(photo_rates[:,1]))
        self.Gamma_HeII = interp.interp1d(redshifts, np.log10(photo_rates[:,2]))
        self.Eps_HI = interp.interp1d(redshifts, np.log10(photo_heat[:,0]))
        self.Eps_HeI = interp.interp1d(redshifts, np.log10(photo_heat[:,1]))
        self.Eps_HeII = interp.interp1d(redshifts, np.log10(photo_heat[:,2]))

    def gHe0(self,redshift):
        """Get photo rate for neutral Helium"""
        log1z = np.log10(1+redshift)
        return 10**self.Gamma_HeI(log1z)

    def gHep(self,redshift):
        """Get photo rate for singly ionized Helium"""
        log1z = np.log10(1+redshift)
        return 10**self.Gamma_HeII(log1z)

    def gH0(self,redshift):
        """Get photo rate for neutral Hydrogen"""
        log1z = np.log10(1+redshift)
        return 10**self.Gamma_HI(log1z)

    def epsHe0(self,redshift):
        """Get photo heating rate for neutral Helium"""
        log1z = np.log10(1+redshift)
        return 10**self.Eps_HeI(log1z)

    def epsHep(self,redshift):
        """Get photo heating rate for singly ionized Helium"""
        log1z = np.log10(1+redshift)
        return 10**self.Eps_HeII(log1z)

    def epsH0(self,redshift):
        """Get photo heating rate for neutral Hydrogen"""
        log1z = np.log10(1+redshift)
        return 10**self.Eps_HI(log1z)

class CoolingRatesKWH92(object):
    """The cooling rates from KWH92, in erg s^-1 cm^-3 (cgs).
    All rates are divided by the abundance of the ions involved in the interaction.
    So we are computing the cooling rate divided by n_e n_X. Temperatures in K.
    None of these rates are original to KWH92, but are taken from Cen 1992,
    and originally from older references. The hydrogen rates in particular are probably inaccurate.
    Cen 1992 modified (arbitrarily) the excitation and ionisation rates for high temperatures.
    There is no collisional excitation rate for He0 - not sure why.
    References:
        Black 1981, from Lotz 1967, Seaton 1959, Burgess & Seaton 1960.
        Recombination rates are from Spitzer 1978.
        Free-free: Spitzer 1978.
    Collisional excitation and ionisation cooling rates are merged.
    """
    def __init__(self, tcmb=2.7255, t5_corr=1e5, recomb=None):
        self.tcmb = tcmb
        if recomb is None:
            self.recomb = RecombRatesCen92()
        else:
            self.recomb = recomb()
        self.t5_corr = t5_corr
        #1 eV in ergs
        self.eVinergs = 1.60218e-12
        #boltzmann constant in erg/K
        self.kB = 1.38064852e-16

    def _t5(self, temp):
        """Commonly used Cen 1992 correction factor for large temperatures.
        This is implemented so that the cooling rates have the right
        asymptotic behaviour. However, Cen erroneously imposes this correction at T=1e5,
        which is too small: the Black 1981 rates these are based on should be good
        until 5e5 at least, where the correction factor has a 10% effect already.
        More modern tables thus impose it at T=5e7, which is still arbitrary but should be harmless.
        """
        return 1+(temp/self.t5_corr)**0.5

    def CollisionalExciteH0(self, temp):
        """Collisional excitation cooling rate for n_H0 and n_e. Gadget calls this BetaH0."""
        return 7.5e-19 * np.exp(-118348.0/temp) /self._t5(temp)

    def CollisionalExciteHeP(self, temp):
        """Collisional excitation cooling rate for n_He+ and n_e. Gadget calls this BetaHep."""
        return 5.54e-17 * temp**(-0.397)*np.exp(-473638./temp)/self._t5(temp)

    def CollisionalExciteHe0(self, temp):
        """This is listed in Cen 92 but neglected in KWH 97, presumably because it is very small."""
        #return 0
        return 9.1e-27 * temp**(-0.1687) * np.exp(-473638/temp) / self._t5(temp)

    def CollisionalIonizeH0(self, temp):
        """Collisional ionisation cooling rate for n_H0 and n_e. Gadget calls this GammaeH0."""
        #Ionisation potential of H0
        return 13.5984 * self.eVinergs * self.recomb.GammaeH0(temp)

    def CollisionalIonizeHe0(self, temp):
        """Collisional ionisation cooling rate for n_H0 and n_e. Gadget calls this GammaeHe0."""
        return 24.5874 * self.eVinergs * self.recomb.GammaeHe0(temp)

    def CollisionalIonizeHeP(self, temp):
        """Collisional ionisation cooling rate for n_H0 and n_e. Gadget calls this GammaeHep."""
        return 54.417760 * self.eVinergs * self.recomb.GammaeHep(temp)

    def CollisionalH0(self, temp):
        """Total collisional cooling for H0"""
        return self.CollisionalExciteH0(temp) + self.CollisionalIonizeH0(temp)

    def CollisionalHe0(self, temp):
        """Total collisional cooling for H0"""
        return self.CollisionalExciteHe0(temp) + self.CollisionalIonizeHe0(temp)

    def CollisionalHeP(self, temp):
        """Total collisional cooling for H0"""
        return self.CollisionalExciteHeP(temp) + self.CollisionalIonizeHeP(temp)

    def RecombHp(self, temp):
        """Recombination cooling rate for H+ and e. Gadget calls this AlphaHp."""
        return 0.75 * self.kB * temp * self.recomb.alphaHp(temp)

    def RecombHeP(self, temp):
        """Recombination cooling rate for He+ and e. Gadget calls this AlphaHep."""
        #I'm not sure why they use 0.75 kT as the free energy of an electron.
        #I would guess this is explained in Spitzer 1978.
        return 0.75 * self.kB * temp * self.recomb.alphaHep(temp)+ self._RecombDielect(temp)

    def RecombHePP(self, temp):
        """Recombination cooling rate for He++ and e. Gadget calls this AlphaHepp."""
        return 0.75 * self.kB * temp * self.recomb.alphaHepp(temp)

    def _RecombDielect(self, temp):
        """Dielectric recombination rate for He+ and e. Gadget calls this Alphad."""
        #What is this magic number?
        return 6.526e-11*self.recomb.alphad(temp)

    def FreeFree(self, temp, zz):
        """Free-free cooling rate for electrons scattering on ions without being captured.
        Factors here are n_e and total ionized species:
            (FreeFree(zz=1)*(n_H+ + n_He+) + FreeFree(zz=2)*n_He++)"""
        return 1.426e-27*np.sqrt(temp)*zz**2*self._gff(temp,zz)

    def _gff(self, temp, zz):
        """Formula for the Gaunt factor. KWH takes this from Spitzer 1978."""
        _ = zz
        return 1.1+0.34*np.exp(-(5.5 - np.log10(temp))**2/3.)

    def InverseCompton(self, temp, redshift):
        """Cooling rate for inverse Compton from the microwave background.
        Multiply this only by n_e. Note the CMB temperature is hardcoded in KWH92 to 2.7."""
        tcmb_red = self.tcmb * (1+redshift)
        #Thompson cross-section in cm^2
        sigmat = 6.6524e-25
        #Radiation density constant, 4 sigma_stefan-boltzmann / c in erg cm^-3 K^-4
        rad_dens = 7.5657e-15
        #Electron mass in g
        me = 9.10938e-28
        #Speed of light in cm/s
        cc = 2.99792e10
        return 4 * sigmat * rad_dens / (me*cc) * tcmb_red**4 * self.kB * (temp - tcmb_red)

class CoolingRatesSherwood(CoolingRatesKWH92):
    """The cooling rates used in the Sherwood simulation, Bolton et al 2017, in erg s^-1 cm^-3 (cgs).
    Differences from KWH92 are updated recombination and collisional ionization rates, and the use of a
    larger temperature correction factor than Cen 92.
    """
    def __init__(self, tcmb=2.7255, recomb=None):
        CoolingRatesKWH92.__init__(self, tcmb = tcmb, t5_corr = 5e7, recomb=RecombRatesVerner96)

class CoolingRatesNyx(CoolingRatesKWH92):
    """The cooling rates used in the Nyx paper Lukic 2014, 1406.6361, in erg s^-1 cm^-3 (cgs).
    All rates are divided by the abundance of the ions involved in the interaction.
    So we are computing the cooling rate divided by n_e n_X. Temperatures in K.
    Major differences from KWH are the use of the Scholz & Walter 1991
    hydrogen collisional cooling rates, a less aggressive high temperature correction for helium, and
    Shapiro & Kang 1987 for free free.
    Older Black 1981 recombination cooling rates are used!
    They use the recombination rates from Verner & Ferland 96, but do not change the cooling rates to match.
    Ditto the ionization rates from Voronov 1997: they should also use these rates for collisional ionisation,
    although this is harder because Sholz & Walter don't break their rates into ionization and excitation.
    References:
        Scholz & Walters 1991 (0.45% accuracy)
        Black 1981 (recombination and helium)
        Shapiro & Kang 1987
    """
    def __init__(self, tcmb=2.7255, recomb=None):
        CoolingRatesKWH92.__init__(self, tcmb = tcmb, t5_corr = 5e7, recomb=recomb)

    def CollisionalH0(self, temp):
        """Collisional cooling rate for n_H0 and n_e. Gadget calls this BetaH0 + GammaeH0.
        Formula from Eq. 23, Table 4 of Scholz & Walters, claimed good to 0.45 %.
        Note though that they have two datasets which differ by a factor of two.
        Differs from Cen 92 by a factor of two."""
        #Technically only good for T > 2000.
        y = np.log(temp)
        #Constant is 0.75/k_B in Rydberg
        Ryd = 2.1798741e-11
        tot = -0.75/self.kB*Ryd/temp
        coeffslowT = [213.7913, 113.9492, 25.06062, 2.762755, 0.1515352, 3.290382e-3]
        coeffshighT = [271.25446, 98.019455, 14.00728, 0.9780842, 3.356289e-2, 4.553323e-4]
        for j in range(6):
            tot += ((temp < 1e5)*coeffslowT[j]+(temp >=1e5)*coeffshighT[j])*(-y)**j
        return 1e-20 * np.exp(tot)

    def RecombHp(self, temp):
        """Recombination cooling rate for H+ and e. Gadget calls this AlphaHp.
        Differs by O(10%) until 3x10^6."""
        return 2.851e-27 * np.sqrt(temp) * (5.914 - 0.5 * np.log(temp) + 0.01184 * temp**(1./3))

    def RecombHePP(self, temp):
        """Recombination cooling rate for H+ and e. Gadget calls this AlphaHepp.
        Differs from Cen 92 by 10% until ~10^7"""
        return 1.140e-26 * np.sqrt(temp) * (6.607 - 0.5 * np.log(temp) + 7.459e-3 * temp**(1./3))

    def _gff(self, temp, zz):
        """Formula for the Gaunt factor from Shapiro & Kang 1987. ZZ is 1 for H+ and He+ and 2 for He++.
        This is almost identical to the KWH rate but not continuous."""
        #This is not continuous. Check the original reference.
        little = (temp/zz**2 <= 3.2e5)
        lt = np.log10(temp/zz**2)
        return little * (0.79464 + 0.1243*lt) + np.logical_not(little) * ( 2.13164 - 0.1240 * lt)

#Fixed point optimization routines from scipy, modified to enforce positivity and use absolute error
from scipy._lib._util import _asarray_validated, _lazywhere

def _del2(p0, p1, d):
    """del2 convergence accelerator"""
    pp = p0 - np.square(p1 - p0) / d
    return np.maximum(pp, np.zeros_like(pp))

def _fixed_point_helper(func, x0, args, xtol, maxiter, use_accel):
    """Helper function from scipy optimize"""
    p0 = x0
    for i in range(maxiter):
        p1 = func(p0, *args)
        if np.all(np.abs(p1-p0) < xtol/2.):
            return p1
        if use_accel and i < maxiter//2:
            p2 = func(p1, *args)
            d = p2 - 2.0 * p1 + p0
            p = _lazywhere(d != 0, (p0, p1, d), f=_del2, fillvalue=p2)
        else:
            p = p1
        p0 = p
    msg = "Failed to converge after %d iterations, value is %s" % (maxiter, p)
    raise RuntimeError(msg)


def fixed_point(func, x0, args=(), xtol=1e-8, maxiter=500, method='del2'):
    """
    Find a fixed point of the function.

    Given a function of one or more variables and a starting point, find a
    fixed-point of the function: i.e. where ``func(x0) == x0``.

    Parameters
    ----------
    func : function
        Function to evaluate.
    x0 : array_like
        Fixed point of function.
    args : tuple, optional
        Extra arguments to `func`.
    xtol : float, optional
        Convergence tolerance, defaults to 1e-08.
    maxiter : int, optional
        Maximum number of iterations, defaults to 500.
    method : {"del2", "iteration"}, optional
        Method of finding the fixed-point, defaults to "del2"
        which uses Steffensen's Method with Aitken's ``Del^2``
        convergence acceleration [1]_. The "iteration" method simply iterates
        the function until convergence is detected, without attempting to
        accelerate the convergence.

    References
    ----------
    .. [1] Burden, Faires, "Numerical Analysis", 5th edition, pg. 80

    Examples
    --------
    >>> from scipy import optimize
    >>> def func(x, c1, c2):
    ...    return np.sqrt(c1/(x+c2))
    >>> c1 = np.array([10,12.])
    >>> c2 = np.array([3, 5.])
    >>> optimize.fixed_point(func, [1.2, 1.3], args=(c1,c2))
    array([ 1.4920333 ,  1.37228132])

    """
    use_accel = {'del2': True, 'iteration': False}[method]
    x0 = _asarray_validated(x0, as_inexact=True)
    return _fixed_point_helper(func, x0, args, xtol, maxiter, use_accel)
