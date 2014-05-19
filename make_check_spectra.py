# -*- coding: utf-8 -*-
"""Generate some velocity widths through DLAs for the checking script"""
import gridspectra as gs
import randspectra as rs
import spectra as ss
import sys
import os.path as path
import numpy as np
import halospectra as hs

np.seterr(all='raise')
np.seterr(under='warn')


def make_stuff(halo):
    """Get the various arrays we want and save them"""
    halo.get_col_density("H",1)
    halo.get_observer_tau("Si",2, force_recompute=True)
    #SiII 1260
    halo.get_tau("Si",2,1260, force_recompute=True)
    halo.get_tau("Si",2,1526, force_recompute=True)
    halo.get_tau("H",1,1215, force_recompute=True)
    halo.get_col_density("Si",2, force_recompute=True)
    halo.get_col_density("Z",-1)
    halo.get_col_density("H",-1)
    halo.save_file()

snapnum=3
sim=7

#Box
base=path.expanduser("~/data/Cosmo/Cosmo5_V6/L10n512/output")
halo = gs.GridSpectra(snapnum, base, numlos=5000)

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

halo = hs.HaloSpectra(snapnum, base, minpart=3000, savefile="halo_spectra_2.hdf5")
make_stuff(halo)

#Tescari self-shielding condition:
#add to _read_particle_data in spectra.py:
#        #Special case H1:
#        if elem == 'H' and ion == 1:
#            # Neutral hydrogen mass frac
#            elem_den *= star.get_reproc_HI(data)[ind]
#            ind3 = np.where(den > 0.1)
#            elem_den = elem_den[ind3]
#            pos = pos[ind3]
#            hh = hh[ind3]
#            if get_tau:
#                temp = temp[ind3]
#                vel = vel[ind3]
#        elif ion != -1:
#            #Cloudy density in physical H atoms / cm^3
#            ind2 = np.where((elem_den > 0)*(den > 0.1))
# halo = gs.GridSpectra(snapnum, base, numlos=1000, savefile="grid_spectra_DLA_noshield.hdf5")
# make_stuff(halo)

# SiII fraction given by CLOUDY from the metallicity and column density.
class ColdenSpectra(ss.Spectra):
    """Spectra with the SiII fraction given by the metallicity and column density."""
    def _get_elem_den(self, elem, ion, den, temp, data, ind, ind2, star):
        """Get the density in an elemental species."""
        zzz = self.get_mass_frac("Z", data, ind)[ind2]/0.76*self.solar[elem]/self.solarz
        return zzz*star.get_reproc_HI(data)[ind][ind2]*np.float32(self.cloudy_table.ion(elem, ion, den, temp)/self.cloudy_table.ion("H", 1, den, temp))

halo = ColdenSpectra(snapnum, base,None, None, savefile="si_colden_spectra.hdf5")
make_stuff(halo)

class SiHISpectra(ss.Spectra):
    """Spectra with the SiII fraction given by n(SiII)/n(Si) = n(HI)/n(H)."""
    def _get_elem_den(self, elem, ion, den, temp, data, ind, ind2, star):
        """Get the density in an elemental species."""
        zzz = self.get_mass_frac("Z", data, ind)[ind2]*self.solar[elem]/self.solarz
        return star.get_reproc_HI(data)[ind][ind2]

halo = SiHISpectra(snapnum, base,None, None, savefile="SiHI_spectra.hdf5")
make_stuff(halo)
