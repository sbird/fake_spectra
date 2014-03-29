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
    halo.get_observer_tau("Si",2)
    #SiII 1260
    halo.get_tau("Si",2,1260)
    halo.get_tau("Si",2,1526)
    halo.get_tau("H",1,1215)
    halo.get_col_density("Si",2)
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

halo = hs.HaloSpectra(snapnum, base, minpart=3000, savefile="halo_spectra.hdf5")
make_stuff(halo)
