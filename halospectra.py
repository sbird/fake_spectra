"""Class to gather and analyse various metal line statistics"""

import numpy as np
import hdfsim
import halocat
import spectra

class HaloSpectra(spectra.Spectra):
    """Generate metal line spectra from simulation snapshot"""
    def __init__(self,num, base, minpart = 400, nbins = 1024, cloudy_dir="/home/spb/codes/ArepoCoolingTables/tmp_spb/"):
        #Load halo centers to push lines through them
        f = hdfsim.get_file(num, base, 0)
        self.OmegaM = f["Header"].attrs["Omega0"]
        self.box = f["Header"].attrs["BoxSize"]
        self.npart=f["Header"].attrs["NumPart_Total"]+2**32*f["Header"].attrs["NumPart_Total_HighWord"]
        min_mass = self.min_halo_mass(minpart)
        f.close()
        (ind, self.sub_mass, cofm, self.sub_radii) = halocat.find_wanted_halos(num, base, min_mass)
        self.NumLos = np.size(self.sub_mass)
        #Random integers from [1,2,3]
        axis = np.random.random_integers(3, size = self.NumLos)
        spectra.Spectra.__init__(self,num, base, cofm, axis, nbins, cloudy_dir)

    def min_halo_mass(self, minpart = 400):
        """Min resolved halo mass in internal Gadget units (1e10 M_sun)"""
        #This is rho_c in units of h^-1 1e10 M_sun (kpc/h)^-3
        rhom = 2.78e+11* self.OmegaM / 1e10 / (1e3**3)
        #Mass of an SPH particle, in units of 1e10 M_sun, x omega_m/ omega_b.
        target_mass = self.box**3 * rhom / self.npart[0]
        min_mass = target_mass * minpart
        return min_mass

