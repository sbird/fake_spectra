# -*- coding: utf-8 -*-
import halospectra as hs
import randspectra as rs
import sys

snapnum=sys.argv[1]
sim=sys.argv[2]
#base="/n/hernquistfs1/mvogelsberger/projects/GFM/Production/Cosmo/Cosmo"+str(sim)+"_V6/L25n512/output/"
#savedir="/n/home11/spb/scratch/Cosmo/Cosmo"+str(sim)+"_V6_512/snapdir_"+str(snapnum).rjust(3,'0')
base="/home/spb/data/Cosmo/Cosmo"+str(sim)+"_V6/L25n256"
savedir="/home/spb/scratch/Cosmo/Cosmo"+str(sim)+"_V6/snapdir_"+str(snapnum).rjust(3,'0')
#halo = hs.HaloSpectra(snapnum, base,3, savefile="halo_spectra_DLA.hdf5", savedir=savedir)
halo = rs.RandSpectra(snapnum, base,numlos=3000,savedir=savedir, savefile="rand_spectra_DLA.hdf5")
halo.get_tau("Si",2,2)
halo.get_tau("H",1,1)
halo.get_col_density("Z",-1)
halo.get_col_density("H",-1)
halo.save_file()

