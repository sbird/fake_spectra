"""Unit system for the spectral code."""
import math
import numpy as np

class UnitSystem(object):
    """Class to store the various physical constants and units that are relevant here. Factored out of Spectra."""
    def __init__(self, UnitMass_in_g=1.98892e43, UnitLength_in_cm=3.085678e21, UnitVelocity_in_cm_per_s = 1e5):
        #Internal gadget mass unit: 1e10 M_sun/h in g/h
        self.UnitMass_in_g = UnitMass_in_g
        #Internal gadget length unit: 1 kpc/h in cm/h
        self.UnitLength_in_cm = UnitLength_in_cm
        #Some constants and unit systems
        self.UnitDensity_in_cgs = self.UnitMass_in_g/self.UnitLength_in_cm**3
        #Internal velocity unit : 1 km/s in cm/s
        self.UnitVelocity_in_cm_per_s = UnitVelocity_in_cm_per_s
        #Internal energy is in erg/g = 1 (km/s)**2 in (cm/s)**2
        self.UnitInternalEnergy_in_cgs = self.UnitVelocity_in_cm_per_s**2
        #Speed of light in cm/s
        self.light = 2.99e10
        #proton mass in g
        self.protonmass=1.67262178e-24
        #Boltzmann constant (cgs)
        self.boltzmann=1.38066e-16
        #Newton's constant in cm^3/g/s^2
        self.gravcgs = 6.674e-8
        #100 km/s/Mpc in 1/s
        self.h100=3.2407789e-18
        #Gas equation of state
        self.gamma=5./3

    def absorption_distance(self, speclen, red):
        """
        Compute X(z), the absorption distance per sightline (dimensionless)
        X(z) = int (1+z)^2 H_0 / H(z) dz
        When dz is small, dz ~ H(z)/c dL, so
        X(z) ~ (1+z)^2 H_0/c dL
        Arguments:
            speclen - spectral length (usually box size in comoving kpc/h)
            red - redshift
        """
        #Units: h/s   s/cm                 kpc/h      cm/kpc
        return self.h100/self.light*speclen*self.UnitLength_in_cm*(1+red)**2

    def redshift_distance(self, speclen, red,omegam0):
        """Compute dz over the box, dz = H(z)/c dL
        Arguments:
            speclen - spectral length (usually box size in comoving kpc/h)
            red - redshift
        """
        #Units: h/s   s/cm                 kpc/h      cm/kpc
        return self.hubble(red, omegam0)/self.light*speclen*self.UnitLength_in_cm

    def hubble(self, z, omegam0):
        """Hubble parameter"""
        return self.h100*np.sqrt(omegam0*(1+z)**3 + (1-omegam0))

    def rho_crit(self, hubble):
        """Get the critical density at z=0 in units of g cm^-3"""
        #H in units of 1/s
        h100=self.h100*hubble
        rho_crit=3*h100**2/(8*math.pi*self.gravcgs)
        return rho_crit
