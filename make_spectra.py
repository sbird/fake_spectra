# -*- coding: utf-8 -*-
"""Generate some velocity widths through DLAs"""
import gridspectra as gs
import randspectra as rs
import spectra as ss
import sys
import os.path as path
import numpy as np

np.seterr(all='raise')
np.seterr(under='warn')


snapnum=sys.argv[1]
sim=sys.argv[2]
#base="/n/hernquistfs1/mvogelsberger/projects/GFM/Production/Cosmo/Cosmo"+str(sim)+"_V6/L25n512/output/"
#savedir="/n/home11/spb/scratch/Cosmo/Cosmo"+str(sim)+"_V6_512/snapdir_"+str(snapnum).rjust(3,'0')
base=path.expanduser("~/data/Cosmo/Cosmo"+str(sim)+"_V6/L25n512/output")
if len(sys.argv) > 3:
    halo = rs.RandSpectra(snapnum, base, numlos=5000, thresh=0)
else:
    halo = gs.GridSpectra(snapnum, base, numlos=5000)

# halo = ss.Spectra(snapnum, base, None, None,savefile = "grid_spectra_DLA.hdf5")
halo.get_col_density("H",1)
halo.get_observer_tau("Si",2)
#SiII 1260
halo.get_tau("Si",2,1260)
halo.get_tau("Si",2,1526)
halo.get_tau("H",1,1215)
halo.get_col_density("Si",2)
halo.get_col_density("Z",-1)
halo.get_col_density("H",-1)
halo.find_nearby_halos()
halo.get_velocity("H",1)
halo.save_file()

