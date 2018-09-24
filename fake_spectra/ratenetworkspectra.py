"""Modified versions of gas properties and spectra that use the rate network."""

import numpy as np
from scipy.interpolate import interp2d
from fake_spectra import gas_properties
from fake_spectra import spectra
from rate_network import RateNetwork

def RateNetworkGasTest():
    """Test that the spline is working."""
    gasprop = RateNetworkGas(3, None)
    dlim = (np.log(1e-7), np.log(3))
    elim = (np.log(20), np.log(3e6))
    randd = (dlim[1] - dlim[0]) * np.random.random(size=2000) + dlim[0]
    randi = (elim[1] - elim[0]) * np.random.random(size=2000) + elim[0]
    spline = gasprop.build_interp(dlim, elim)
    for dd, ii in zip(randd, randi):
        spl = spline(dd, ii)[0]
        rate = np.log(gasprop.rates.get_neutral_fraction(np.exp(dd), np.exp(ii)))
        assert np.abs(spl - rate) < 1e-5

class RateNetworkGas(gas_properties.GasProperties):
    """Replace the get_reproc_HI function with something that solves the rate network. Optionally can also do self-shielding."""
    def __init__(self, redshift, absnap, hubble=0.71, fbar=0.17, units=None, sf_neutral=True, selfshield=False, photo_factor=1):
        super().__init__(redshift, absnap, hubble=hubble, fbar=fbar, units=units, sf_neutral=sf_neutral)
        self.rates = RateNetwork(redshift, photo_factor = photo_factor, f_bar = fbar, cool="KWH", recomb="C92", selfshield=selfshield, treecool_file="TREECOOL_hm_2012_sherwood")
        self.temp_factor = 1
        self.gamma_factor = 1

    def build_interp(self, dlim, elim, sz=500):
        """Build the interpolator"""
        #Build interpolation
        densgrid = np.linspace(dlim[0], dlim[1], 2*sz)
        ienergygrid = np.linspace(elim[0], elim[1], sz)
        dgrid, egrid = np.meshgrid(densgrid, ienergygrid)
        lh0grid = np.log(self.rates.get_neutral_fraction(np.exp(dgrid), np.exp(egrid)))
        #We assume primordial helium
        spline = interp2d(densgrid, ienergygrid, lh0grid, kind='cubic')
        return spline

    def get_reproc_HI(self, part_type, segment):
        """Get a neutral hydrogen fraction using a rate network which reads temperature and density of the gas."""
        #expecting units of atoms/cm^3
        density = np.log(self.get_code_rhoH(part_type, segment))
        #expecting units of 10^-10 ergs/g
        ienergy = np.log(self.absnap.get_data(part_type, "InternalEnergy", segment=segment)*self.units.UnitInternalEnergy_in_cgs/1e10)
        spline = self.build_interp(dlim=(np.min(density), np.max(density)), elim=(np.min(ienergy), np.max(ienergy)))
        #spline(density, ienergy) evaluates on a Ndens x Nenerg grid!
        nh0 = np.exp([spline(dd, ii)[0] for (dd, ii) in zip(density, ienergy)])
        return nh0

    def _get_ienergy_rescaled(self, density, ienergy, density0):
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
        ienergy *= self.temp_factor
        #Adjust slope by same factor: note use gamma_factor -1 so gamma_factor = 1 means no change.
        if self.gamma_factor != 1.:
            ienergy *= (density/density0)**self.gamma_factor-1.
        return ienergy

class RateNetworkSpectra(spectra.Spectra):
    """Generate spectra with a neutral fraction from a rate network"""
    def __init__(self, *args, photo_factor = 1, selfshield=False, **kwargs):
        super().__init__(*args, **kwargs)
        try:
            self.gasprop = RateNetworkGas(redshift = self.red, absnap = self.snapshot_set, hubble=self.hubble, units = self.units, sf_neutral=False, photo_factor = photo_factor, selfshield=selfshield)
        except AttributeError:
            pass
