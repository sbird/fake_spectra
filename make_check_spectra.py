# -*- coding: utf-8 -*-
"""Generate some velocity widths through DLAs for the checking script"""
import gridspectra as gs
import spectra as ss
import os.path as path
import numpy as np
import halospectra as hs

np.seterr(all='raise')
np.seterr(under='warn')


def make_stuff(halo):
    """Get the various arrays we want and save them"""
    halo.get_density("H",1)
    halo.get_observer_tau("Si",2, force_recompute=True)
    #SiII 1260
    halo.get_tau("Si",2,1260, force_recompute=True)
    halo.get_tau("Si",2,1526, force_recompute=True)
    halo.get_tau("H",1,1215, force_recompute=True)
    halo.get_density("Si",2, force_recompute=True)
    halo.get_density("Z",-1)
    halo.get_density("H",-1)
    halo.save_file()

snapnum=3
sim=7

#Box
base10=path.expanduser("~/data/Cosmo/Cosmo5_V6/L10n512/output")
halo = gs.GridSpectra(snapnum, base10, numlos=5000)

make_stuff(halo)

#Spectral resolution
base=path.expanduser("~/data/Cosmo/Cosmo"+str(sim)+"_V6/L25n512/output")

halo = gs.GridSpectra(snapnum, base, numlos=5000, res=0.5, savefile="grid_spectra_DLA_res.hdf5")
make_stuff(halo)

#Pecvel
#Needs a recompile with peculiar velocities off
# halo = gs.GridSpectra(snapnum, base, numlos=5000, savefile="grid_spectra_DLA_pecvel.hdf5")
# make_stuff(halo)

#Tophat
#Needs a recompile with -DTOP_HAT_KERNEL
# halo = gs.GridSpectra(snapnum, base, numlos=5000, savefile="grid_spectra_DLA_tophat.hdf5")
# make_stuff(halo)

#Attentuation
halo = gs.GridSpectra(snapnum, base, numlos=5000, cdir=path.expanduser("~/codes/cloudy_tables/ion_out_no_atten/"), savefile="grid_spectra_DLA_no_atten.hdf5")
make_stuff(halo)

#Tescari halos

class TescariSpectra(hs.HaloSpectra):
    """Spectra with the SiII fraction given by n(SiII)/n(Si) = n(HI)/n(H)."""
    def _get_elem_den(self, elem, ion, den, temp, data, ind, ind2, star):
        """Get the density in an elemental species. Broken out so it can be over-ridden by child classes."""
        #Make sure temperature doesn't overflow the cloudy table
        if np.max(temp) > 10**8.6:
            temp2 = np.array(temp)
            temp2[np.where(temp2 > 10**8.6)] = 10**8.6
        else:
            temp2 = temp
        ions = np.ones_like(den)
        ind3 = np.where(den < 0.1)
        ions[ind3] = np.float32(self.cloudy_table.ion(elem, ion, den[ind3], temp2[ind3]))
        return ions

#Use a 10 Mpc box for better comparison
halo = TescariSpectra(snapnum, base10, minpart=3000, savefile="halo_spectra_2.hdf5",cdir=path.expanduser("~/codes/cloudy_tables/ion_out_no_atten/"))
make_stuff(halo)

# SiII fraction given by CLOUDY from the metallicity and column density.
class ColdenSpectra(ss.Spectra):
    """Spectra with the SiII fraction given by the metallicity and column density."""
    def _get_elem_den(self, elem, ion, den, temp, data, ind, ind2, star):
        """Get the density in an elemental species."""
        #Make sure temperature doesn't overflow cloudy
        if np.max(temp) > 10**8.6:
            temp2 = np.array(temp)
            temp2[np.where(temp2 > 10**8.6)] = 10**8.6
        else:
            temp2 = temp
        return star.get_reproc_HI(data)[ind][ind2]*np.float32(self.cloudy_table.ion(elem, ion, den, temp2)/self.cloudy_table.ion("H", 1, den, temp2))

    def get_mass_frac(self, elem, data, ind):
        """Get the mass fraction in an elemental species."""
        zzz = ss.Spectra.get_mass_frac(self, "Z", data, ind)/self.solarz*self.solar[elem]/0.76
        return zzz


halo = ColdenSpectra(snapnum, base,None, None, savefile="si_colden_spectra.hdf5")
make_stuff(halo)

class SiHISpectra(ss.Spectra):
    """Spectra with the SiII fraction given by n(SiII)/n(Si) = n(HI)/n(H)."""
    def _get_elem_den(self, elem, ion, den, temp, data, ind, ind2, star):
        """Get the density in an elemental species."""
        return star.get_reproc_HI(data)[ind][ind2]

halo = SiHISpectra(snapnum, base,None, None, savefile="SiHI_spectra.hdf5")
make_stuff(halo)
