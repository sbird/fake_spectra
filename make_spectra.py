# -*- coding: utf-8 -*-
"""Generate some velocity widths through DLAs"""
import randspectra as rs
import sys
import os.path as path

snapnum=sys.argv[1]
sim=sys.argv[2]
#base="/n/hernquistfs1/mvogelsberger/projects/GFM/Production/Cosmo/Cosmo"+str(sim)+"_V6/L25n512/output/"
#savedir="/n/home11/spb/scratch/Cosmo/Cosmo"+str(sim)+"_V6_512/snapdir_"+str(snapnum).rjust(3,'0')
base=path.expanduser("~/data/Cosmo/Cosmo"+str(sim)+"_V6/L25n512")
halo = rs.RandSpectra(snapnum, base)
#halo.save_file()
halo.get_observer_tau("Si",2, force_recompute=True)
halo.get_col_density("H",1)
#halo.get_tau("H",1,1)
halo.get_col_density("Z",-1)
halo.get_col_density("H",-1)
halo.save_file()

