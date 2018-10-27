# -*- coding: utf-8 -*-
"""A rate network for neutral hydrogen following
Katz, Weinberg & Hernquist 1996, eq. 28-32."""
import numpy as np
import scipy.interpolate as interp
import scipy.optimize



class RateNetwork(object):
    """A rate network for neutral hydrogen following
    Katz, Weinberg & Hernquist 1996, astro-ph/9509107, eq. 28-32."""
    def __init__(self,redshift, photo_factor = 1., f_bar = 0.17, converge = 1e-7, selfshield=True, cool="Nyx", recomb="B06", collisional=True, treecool_file="TREECOOL_ep_2018p"):
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
        elif cool == "Nyx":
            self.cool = CoolingRatesNyx()
        else:
            self.cool = CoolingRatesNewGadget()
        #proton mass in g
        self.protonmass = 1.67262178e-24
        self.redshift = redshift
        self.converge = converge
        self.selfshield = selfshield
        self.collisional = collisional
        zz = [0, 1, 2, 3, 4, 5, 6, 7,8]
        #Tables for the self-shielding correction. Note these are not well-measured for z > 5!
        gray_opac = [2.59e-18,2.37e-18,2.27e-18, 2.15e-18, 2.02e-18, 1.94e-18, 1.82e-18, 1.71e-18, 1.60e-18]
        self.Gray_ss = interp.InterpolatedUnivariateSpline(zz, gray_opac)

    def get_temp(self, density, ienergy, helium=0.24):
        """Get the equilibrium temperature at given internal energy.
        density is gas density in protons/cm^3
        Internal energy is in J/kg == 10^-10 ergs/g.
        helium is a mass fraction"""
        ne = self.get_equilib_ne(density, ienergy, helium)
        nh = density * (1-helium)
        return self._get_temp(ne/nh, ienergy, helium)

    def get_cooling_rate(self, density, ienergy, helium=0.24):
        """Get the total cooling rate for a temperature and density."""
        ne = self.get_equilib_ne(density, ienergy, helium)
        nh = density * (1-helium)
        temp = self._get_temp(ne/nh, ienergy, helium)
        nH0 = self._nH0(nh, temp, ne)
        nHe0 = self._nHe0(nh, temp, ne)
        nHep = self._nHep(nh, temp, ne)
        nHp = self._nHp(nh, temp, ne)
        nHep = self._nHep(nh, temp, ne)
        nHepp = self._nHepp(nh, temp, ne)
        #This is the collisional excitation and ionisation rate.
        LambdaCollis = ne * (self.cool.CollisionalH0(temp) * nH0 +
                             self.cool.CollisionalHe0(temp) * nHe0 +
                             self.cool.CollisionalHeP(temp) * nHep)
        LambdaRecomb = ne * (self.cool.RecombHp(temp) * nHp +
                             self.cool.RecombHeP(temp) * nHep +
                             self.cool.RecombHePP(temp) * nHepp)
        LambdaFF = ne * (self.cool.FreeFree(temp, 1)*(nHp + nHep) + self.cool.FreeFree(temp, 2)*nHepp)
        LambdaCmptn = ne * self.cool.InverseCompton(temp, self.redshift)

        return LambdaCollis + LambdaRecomb +  LambdaCmptn + LambdaFF 

    def get_equilib_ne(self, density, ienergy,helium=0.24):
        """Solve the system of equations for photo-ionisation equilibrium,
        starting with ne = nH and continuing until convergence.
        density is gas density in protons/cm^3
        Internal energy is in J/kg == 10^-10 ergs/g.
        helium is a mass fraction.
        """
        #Get hydrogen number density
        nh = density * (1-helium)
        rooted = lambda ne: self._ne(nh, self._get_temp(ne/nh, ienergy, helium=helium), ne, helium=helium)
        ne = scipy.optimize.fixed_point(rooted, nh,xtol=self.converge)
        assert np.all(np.abs(rooted(ne) - ne) < self.converge)
        return ne

    def get_ne_by_nh(self, density, ienergy, helium=0.24):
        """Same as above, but get electrons per proton."""
        return self.get_equilib_ne(density, ienergy, helium)/(density*(1-helium))

    def get_neutral_fraction(self, density, ienergy, helium=0.24):
        """Get the neutral hydrogen fraction at a given temperature and density.
        density is gas density in protons/cm^3
        Internal energy is in J/kg == 10^-10 ergs/g.
        helium is a mass fraction.
        """
        ne = self.get_equilib_ne(density, ienergy, helium=helium)
        nh = density * (1-helium)
        temp = self._get_temp(ne/nh, ienergy, helium)
        return self._nH0(nh, temp, ne) / nh

    def _nH0(self, nh, temp, ne):
        """The neutral hydrogen number density. Eq. 33 of KWH."""
        alphaHp = self.recomb.alphaHp(temp)
        GammaeH0 = self.collisional * self.recomb.GammaeH0(temp)
        photorate = self.photo.gH0(self.redshift)/ne*self.photo_factor*self._self_shield_corr(nh, temp)
        return nh * alphaHp/ (alphaHp + GammaeH0 + photorate)

    def _nHp(self, nh, temp, ne):
        """The ionised hydrogen number density. Eq. 34 of KWH."""
        return nh - self._nH0(nh, temp, ne)

    def _nHep(self, nh, temp, ne):
        """The ionised helium number density, divided by the helium number fraction. Eq. 35 of KWH."""
        alphaHep = self.recomb.alphaHep(temp) + self.recomb.alphad(temp)
        alphaHepp = self.recomb.alphaHepp(temp)
        photofac = self.photo_factor*self._self_shield_corr(nh, temp)
        GammaHe0 = self.collisional * self.recomb.GammaeHe0(temp) + self.photo.gHe0(self.redshift)/ne*photofac
        GammaHep = self.collisional * self.recomb.GammaeHep(temp) + self.photo.gHep(self.redshift)/ne*photofac
        return nh / (1 + alphaHep / GammaHe0 + GammaHep/alphaHepp)

    def _nHe0(self, nh, temp, ne):
        """The neutral helium number density, divided by the helium number fraction. Eq. 36 of KWH."""
        alphaHep = self.recomb.alphaHep(temp) + self.recomb.alphad(temp)
        photofac = self.photo_factor*self._self_shield_corr(nh, temp)
        GammaHe0 = self.collisional * self.recomb.GammaeHe0(temp) + self.photo.gHe0(self.redshift)/ne*photofac
        return self._nHep(nh, temp, ne) * alphaHep / GammaHe0

    def _nHepp(self, nh, temp, ne):
        """The doubly ionised helium number density, divided by the helium number fraction. Eq. 37 of KWH."""
        photofac = self.photo_factor*self._self_shield_corr(nh, temp)
        GammaHep = self.collisional * self.recomb.GammaeHep(temp) + self.photo.gHep(self.redshift)/ne*photofac
        alphaHepp = self.recomb.alphaHepp(temp)
        return self._nHep(nh, temp, ne) * GammaHep / alphaHepp

    def _ne(self, nh, temp, ne, helium=0.24):
        """The electron number density. Eq. 38 of KWH."""
        yy = helium / 4 / (1 - helium)
        return self._nHp(nh, temp, ne) + yy * self._nHep(nh, temp, ne) + 2* yy * self._nHepp(nh, temp, ne)

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

    def _get_temp(self, nebynh, ienergy, helium=0.24):
        """Compute temperature (in K) from internal energy and electron density.
           Uses: internal energy
                 electron abundance per H atom (ne/nH)
                 hydrogen mass fraction (0.76)
           Internal energy is in J/kg, internal gadget units, == 10^-10 ergs/g.
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
        muienergy = 4 / (hy_mass * (3 + 4*nebynh) + 1)*ienergy*1e10
        #Boltzmann constant (cgs)
        boltzmann=1.38066e-16
        gamma=5./3
        #So for T in K, boltzmann in erg/K, internal energy has units of erg/g
        temp = (gamma-1) * self.protonmass / boltzmann * muienergy
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
        Temp in K."""
        return 1.9e-3 / np.power(temp,1.5) * np.exp(-4.7e5/temp)*(1+0.3*np.exp(-9.4e4/temp))

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
    def __init__(self, treecool_file="TREECOOL"):
        #Format of the treecool table:
        # log_10(1+z), Gamma_HI, Gamma_HeI, Gamma_HeII,  Qdot_HI, Qdot_HeI, Qdot_HeII,
        # where 'Gamma' is the photoionization rate and 'Qdot' is the photoheating rate.
        # The Gamma's are in units of s^-1, and the Qdot's are in units of erg s^-1.
        data = np.loadtxt(treecool_file)
        redshifts = data[:,0]
        photo_rates = data[:,1:4]
        assert np.shape(redshifts)[0] == np.shape(photo_rates)[0]
        self.Gamma_HI = interp.InterpolatedUnivariateSpline(redshifts, photo_rates[:,0])
        self.Gamma_HeI = interp.InterpolatedUnivariateSpline(redshifts, photo_rates[:,1])
        self.Gamma_HeII = interp.InterpolatedUnivariateSpline(redshifts, photo_rates[:,2])

    def gHe0(self,redshift):
        """Get photo rate for neutral Helium"""
        log1z = np.log10(1+redshift)
        return self.Gamma_HeI(log1z)

    def gHep(self,redshift):
        """Get photo rate for singly ionized Helium"""
        log1z = np.log10(1+redshift)
        return self.Gamma_HeII(log1z)

    def gH0(self,redshift):
        """Get photo rate for neutral Hydrogen"""
        log1z = np.log10(1+redshift)
        return self.Gamma_HI(log1z)

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
    def __init__(self, tcmb=2.7255, recomb=None):
        self.tcmb = tcmb
        if recomb is None:
            self.recomb = RecombRatesCen92()
        else:
            self.recomb = recomb
        #1 eV in ergs
        self.eVinergs = 1.60218e-12
        #boltzmann constant in erg/K
        self.kB = 1.38064852e-16

    def _t5(self, temp):
        """Commonly used Cen 1992 correction factor for large temperatures."""
        return 1+(temp/1e5)**0.5

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

class CoolingRatesNyx(CoolingRatesKWH92):
    """The cooling rates used in the Nyx paper Lukic 2014, 1406.6361, in erg s^-1 cm^-3 (cgs).
    All rates are divided by the abundance of the ions involved in the interaction.
    So we are computing the cooling rate divided by n_e n_X. Temperatures in K.
    Major differences from KWH are the use of the Scholz & Walter 1991
    hydrogen collisional cooling rates, a less aggressive high temperature correction for helium, and
    Shapiro & Kang 1987 for free free.
    Older Black 1981 recombination cooling rates are used, but I don't know why!
    They use the recombination rates from Voronov 1997, but don't seem to realise that
    this should also change the cooling rates.
    Ditto the ionization rates from Verner & Ferland 96: they should also use these rates for collisional ionisation.
    References:
        Scholz & Walters 1991 (0.45% accuracy)
        Black 1981 (recombination and helium)
        Shapiro & Kang 1987
    """
    def _t5(self, temp):
        """
        Lukic uses a less aggressive correction factor for large temperatures than Cen 1992.
        No explanation is given for this in the paper, but he explained privately that this is probably
        so that it matches Black 1981 for all the regime where that code is valid.
        This factor increases the excitation cooling rate for helium by about a factor of ten and
        thus changes the cooling curve substantially.
        """
        return 1+(temp/5e7)**0.5

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

class CoolingRatesNewGadget(CoolingRatesKWH92):
    """These are the cooling rates that I think at the moment should be used.
    Collisional Ionization and recombination rates are from Voronov 97 and Verner & Ferland 96, respectively.
    Helium excitation is from Cen 1992. Free-free is Spitzer 1978 (as Shapiro & Kang 1987 is discontinuous).
    Hydrogen excitation cooling is Scholz & Walters 1991, which *only* includes the Gamma 1s-2s and
    Gamma 1s-2p terms. However it is not dominant except at high densities where everything is neutral.
    Notably this means that Cen's t5 correction factor only appears in the
    He+ collisional excitation rate, where it should be safely
    negligible (because at high temperatures ionisation dominates).
    """
    def __init__(self, tcmb=2.7255):
        CoolingRatesKWH92.__init__(self, tcmb=tcmb, recomb=RecombRatesVerner96)

    def CollisionalExciteH0(self, temp):
        """Collisional excitation cooling rate for n_H0 and n_e. Gadget calls this BetaH0.
        Formula from Eq. 16, 17, Table 3 of Scholz & Walters 1991.
        This *only* includes the Gamma 1s-2s and Gamma 1s-2p terms. At worst this may underestimate cooling by 20%.
        However it is not dominant except at high densities where everything is neutral, and we probably miss
        molecular lines there anyway.
        """
        y = np.log(temp)
        bblowT = [3.299613e1,1.858848e1, 6.052265, 8.603783e-1, 5.717760e-2, 1.451330e-3]
        cclowT = [1.630155e2,8.795711e1, 2.057117e1, 2.359573, 1.339059e-1,3.021507e-3]
        bbmedT = [2.869759e2, 1.077956e2, 1.524107e1, 1.080538, 3.836975e-2, 5.467273e-4]
        ccmedT = [5.279996e2, 1.939399e2, 2.718982e1, 1.883399, 6.462462e-2, 8.811076e-4]
        bbhighT = [2.7604708e3, 7.9339351e2, 9.1198462e1, 5.1993362, 1.4685343e-1, 1.6404093e-3]
        cchighT = [2.8133632e3, 8.1509685e2, 9.4418414e1, 5.4280565, 1.5467120e-1, 1.7439112e-3]
        gamma2s = 0.
        gamma2p = 0.
        for j in range(6):
            gamma2s += (-1*(temp < 6e4)*bblowT[j]+(temp >=6e4)*(temp < 6e6)*bbmedT[j]-1*(temp > 6e6)*bbhighT[j])*(-y)**j
            gamma2p += (-1*(temp < 6e4)*cclowT[j]+(temp >=6e4)*(temp < 6e6)*ccmedT[j]-1*(temp > 6e6)*cchighT[j])*(-y)**j
        #10.2 eV in erg
        return 10.2 * 1.60184e-12 * (np.exp(gamma2s) + np.exp(gamma2p)) * np.exp(-11606* 10.2/ temp)
