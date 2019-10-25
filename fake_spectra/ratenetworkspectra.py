"""Modified versions of gas properties and spectra that use the rate network."""

import numpy as np
from ._spectra_priv import _interpolate_2d
from . import gas_properties
from . import spectra
from .rate_network import RateNetwork

class RateNetworkGas(gas_properties.GasProperties):
    """Replace the get_reproc_HI function with something that solves the rate network. Optionally can also do self-shielding."""
    def __init__(self, redshift, absnap, hubble=0.71, fbar=0.17, units=None, sf_neutral=True, selfshield=True, photo_factor=1, temp_factor=1, gamma_factor=1):
        super().__init__(redshift, absnap, hubble=hubble, fbar=fbar, units=units, sf_neutral=sf_neutral)
        self.rates = RateNetwork(redshift, photo_factor = photo_factor, f_bar = fbar, selfshield=selfshield, treecool_file="data/TREECOOL_ep_2018p")
        self.temp_factor = temp_factor
        self.gamma_factor = gamma_factor
        self.build_interp(dlim=(-15, 17), elim=(1, 26))

    def build_interp(self, dlim, elim, sz=750):
        """Build the interpolator"""
        #Build interpolation
        self.densgrid = np.linspace(dlim[0], dlim[1], 2*sz)
        self.ienergygrid = np.linspace(elim[0], elim[1], sz)
        dgrid, egrid = np.meshgrid(self.densgrid, self.ienergygrid)
        self.lh0grid = np.log(self.rates.get_neutral_fraction(np.exp(dgrid), np.exp(egrid)))
        #We assume primordial helium

    def get_reproc_HI(self, part_type, segment):
        """Get a neutral hydrogen fraction using a rate network which reads temperature and density of the gas."""
        #expecting units of atoms/cm^3
        density = self.get_code_rhoH(part_type, segment)
        #expecting units of 10^-10 ergs/g
        ienergy = self.absnap.get_data(part_type, "InternalEnergy", segment=segment)*self.units.UnitInternalEnergy_in_cgs/1e10
        ienergy = self._get_ienergy_rescaled(density, ienergy)
        #Correct internal energy to the internal energy of a cold cloud if we are on the star forming equation of state.
        if self.sf_neutral:
            conv = np.float32(self.units.UnitDensity_in_cgs*self.hubble**2/(self.units.protonmass)*(1+self.redshift)**3)
            ind = np.where(density > self.PhysDensThresh/0.76/conv)
            meanweight = 4.0 / (1 + 3 * 0.76)
            EgySpecCold = 1 / (meanweight * (5./3.-1)) * (self.units.boltzmann / self.units.protonmass) * 1000
            ienergy[ind] = EgySpecCold/1e10
        density = np.log(density)
        ienergy = np.log(ienergy)
        if (np.max(self.densgrid) < np.max(density)) or (np.min(self.densgrid) > np.min(density)):
            raise ValueError("Density out of range: interp %g -> %g. Present: %g -> %g" % (np.min(self.densgrid), np.max(self.densgrid), np.min(density), np.max(density)))
        if (np.max(self.ienergygrid) < np.max(ienergy)) or (np.min(self.ienergygrid) > np.min(ienergy)):
            raise ValueError("Ienergy out of range: interp %g -> %g. Present: %g -> %g" % (np.min(self.ienergygrid), np.max(self.ienergygrid), np.min(ienergy), np.max(ienergy)))
        nH0 = np.exp(_interpolate_2d(density, ienergy, self.densgrid, self.ienergygrid, self.lh0grid))
        return nH0

    def _get_ienergy_rescaled(self, density, ienergy):
        """Get the internal energy, rescaled to give the desired equation of state.
        Technically the e. of s. normally used is:
            T = T_0 (rho / rho_0)^(gamma-1)
        However in photoionisation equilibrium the electron density depends very weakly
        on the temperature, and so T/T_0 = U/U_0
        So we can just rescale the internal energy:
        when T_0 -> T_0' U -> U * T_0'/T_0.
        Ditto for gamma, when gamma -> gamma' we have:
        U -> U (rho/rho_0) ^(gamma'-gamma)
        Note this means that if any particle lies off the original equation of state,
        it lies off the new one by a similar amount; the dispersion is preserved!
        """
        #Adjust temperature by desired factor, to give desired equation of state.
        omegab = 0.0445
        rhoc = self.units.rho_crit(self.hubble) * (1+self.redshift)**3
        overden = self.units.protonmass * density /(omegab * rhoc)
        ienergy *= self.temp_factor
        #Adjust slope by same factor: note use gamma_factor -1 so gamma_factor = 1 means no change.
        if self.gamma_factor != 1.:
            ienergy *= (overden)**self.gamma_factor-1.
        return ienergy

class RateNetworkSpectra(spectra.Spectra):
    """Generate spectra with a neutral fraction from a rate network"""
    def __init__(self, *args, photo_factor = 1, sf_neutral=True, selfshield=True, **kwargs):
        kwargs["gasprop"]=RateNetworkGas
        kwargs["sf_neutral"] = sf_neutral
        kwargs["gasprop_args"] = {"photo_factor": photo_factor, "selfshield" : selfshield}
        super().__init__(*args, **kwargs)
