"""Modified versions of gas properties and spectra that use the rate network."""

import numpy as np
from ._spectra_priv import _interpolate_2d
from . import gas_properties
from . import spectra
from .rate_network import RateNetwork

class RateNetworkGas(gas_properties.GasProperties):
    """Replace the get_reproc_HI function with something that solves the rate network. Optionally can also do self-shielding."""
    def __init__(self, redshift, absnap, hubble=0.71, fbar=0.17, units=None, sf_neutral=True, temp_factor=1, gamma_factor=1, **kwargs):
        super().__init__(redshift, absnap, hubble=hubble, fbar=fbar, units=units, sf_neutral=sf_neutral)
        self.rates = RateNetwork(redshift, f_bar = fbar, **kwargs)
        self.temp_factor = temp_factor
        self.gamma_factor = gamma_factor
        self.maxdens = self.PhysDensThresh/0.76
        dmax = 5
        dsz=1000
        if self.sf_neutral:
            dmax = np.log(self.maxdens)
            dsz = 500
        self.build_interp(dlim=(-16, dmax), elim=(2, 21),tsz=500, dsz=dsz)

    def build_interp(self, dlim, elim, tsz=500, dsz=1000):
        """Build the interpolator"""
        #Build interpolation
        self.densgrid = np.linspace(dlim[0], dlim[1], dsz)
        self.ienergygrid = np.linspace(elim[0], elim[1], tsz)
        dgrid, egrid = np.meshgrid(self.densgrid, self.ienergygrid)
        self.lh0grid = np.zeros_like(dgrid)
        #We assume primordial helium
        for i in range(dsz):
            self.lh0grid[:,i] = np.log(self.rates.get_neutral_fraction(np.exp(dgrid[:,i]), np.exp(egrid[:,i])))

    def get_reproc_HI(self, part_type, segment):
        """Get a neutral hydrogen fraction using a rate network which reads temperature and density of the gas."""
        #expecting units of atoms/cm^3
        density = self.get_code_rhoH(part_type, segment)
        #expecting units of 10^-10 ergs/g
        ienergy = self.absnap.get_data(part_type, "InternalEnergy", segment=segment)*self.units.UnitInternalEnergy_in_cgs/1e10
        ienergy = self._get_ienergy_rescaled(density, ienergy)
        ldensity = np.log(density)
        lienergy = np.log(ienergy)
        #Clamp the temperatures : hot gas has the same neutral fraction of 0 anyway.
        ie = np.where(lienergy >= np.max(self.ienergygrid))
        lienergy[ie] = np.max(self.ienergygrid)*0.99
        ie = np.where(lienergy <= np.min(self.ienergygrid))

        lienergy[ie] = np.min(self.ienergygrid)*1.01

        nH0 = np.ones_like(density)
        ii = np.where(ldensity < np.max(self.densgrid))
        if (np.max(self.ienergygrid) < np.max(lienergy[ii])) or (np.min(self.ienergygrid) > np.min(lienergy[ii])):
            raise ValueError("Ienergy out of range: interp %g -> %g. Present: %g -> %g" % (np.min(self.ienergygrid), np.max(self.ienergygrid), np.min(lienergy[ii]), np.max(lienergy[ii])))
        #Correct internal energy to the internal energy of a cold cloud if we are on the star forming equation of state.
        nH0[ii] = np.exp(_interpolate_2d(ldensity[ii], lienergy[ii], self.densgrid, self.ienergygrid, self.lh0grid))
        ii2 = np.where(ldensity >= np.max(self.densgrid))
        if np.size(ii2) > 0:
            if self.sf_neutral:
                if self.redshift_coverage:
                    ssnH0 = self._neutral_fraction(density[ii2], 1e4)
                    nH0[ii2] = ssnH0
                else:
                    nH0[ii2] = 1.
            else:
                nH0[ii2] = self.rates.get_neutral_fraction(density[ii2], ienergy[ii2])

        assert np.all(np.logical_not(np.isnan(nH0)))

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

            ienergy *= (overden)**(self.gamma_factor-1.)
        assert np.all(np.logical_not(np.isnan(ienergy)))
        assert np.all(ienergy > 0)

        return ienergy

class RateNetworkSpectra(spectra.Spectra):
    """Generate spectra with a neutral fraction from a rate network"""
    def __init__(self, *args, photo_factor = 1, sf_neutral=True,
            selfshield=True, temp_factor = 1, gamma_factor = 1,
            hubble = 0.71, fbar = 0.17,
            treecool_file = "data/TREECOOL_ep_2018p", **kwargs):
        kwargs["gasprop"]=RateNetworkGas
        kwargs["sf_neutral"] = sf_neutral
        kwargs["gasprop_args"] = {"photo_factor" : photo_factor,
            "selfshield" : selfshield, "temp_factor" : temp_factor,
            "gamma_factor" : gamma_factor, "hubble" : hubble,
            "fbar" : fbar, "treecool_file" : treecool_file}
        super().__init__(*args, **kwargs)
