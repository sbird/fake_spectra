import fake_spectra
#print("fake_spectra imported")
from fake_spectra.randspectra import RandSpectra
#print("haloassigned_spectra imported")
from fake_spectra.plot_spectra import PlottingSpectra
#print("PlottingSpectra imported")
from mpi4py import MPI
import numpy as np
import time

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

tss = time.asctime()
print('Rank =', rank, 'started at :', tts, flush=True)

numlos = 10000
ndla = 1000
thresh_t = 10**20.3
### the array to store added column density in

rank_str = str(rank)
rr = RandSpectra(34, "/rhome/mqezl001/bigdata/TNG/TNG100-1/output/snapdir_034/" +rank_str + "/", MPI, comm,  thresh = 0.0, kernel='tophat',ndla = 1000, numlos=1000,savedir="/rhome/mqezl001/bigdata/TNG/TNG100-1/postprocessing/randspectra/Snap_034/parallel", savefile="spectra_34."+rank_str+".hdf5")

#### Calculate spectra for 100 hdf5 files
rr.get_tau("H",1,1215)
#Lyman-beta
rr.get_tau("H",1,1025)
rr.get_col_density("H",1)
rr.get_col_density("H",-1)
rr.get_col_density("Z", -1)

rr.save_file()
del rr

tsd = time.asctime()
print('Rank =', rank, 'was Done at :', tsd, flush = True)
