import fake_spectra_parallel
from fake_spectra_parallel.randspectra import RandSpectra
from fake_spectra_parallel.plot_spectra import PlottingSpectra
from fake_spectra_parallel.ratenetworkspectra import RateNetworkGas
from mpi4py import MPI
import numpy as np
import time
import argparse

def get_spec(num, base, savefile, savedir, res, numlos, thresh):
   comm = MPI.COMM_WORLD
   rank = comm.Get_rank()
   size = comm.Get_size()

   tss = time.asctime()
   print('Rank =', rank, 'started!', tss, flush=True)
   
   rr = RandSpectra(num = num, numlos=numlos, ndla=numlos, thresh=thresh, res=res, base= base , MPI = MPI, comm = comm, savedir=savedir, savefile='spectra_'+str(rank)+'.hdf5', gasprop=RateNetworkGas, gasprop_args={"selfshield":True, "treecool_file":"./TREECOOL_ep_2018p", "cool":"KWH","recomb":"Cen92"}, kernel='tophat')

   rr.get_tau("H",1,1215)
   rr.get_col_density("H",1)
   rr.save_file()

   tsd = time.asctime() 
   print('Rank = ', rank, 'is done !', tsd, flush=True)
   del rr

if __name__ == '__main__':

   parser = argparse.ArgumentParser()
   parser.add_argument('-num', type=int, required=True)
   parser.add_argument('-base', type=str, required=True)
   parser.add_argument('-savedir', type=str, required=True)
   parser.add_argument('-res',  type=float, required=True)
   parser.add_argument('-numlos', type=int, required=True)
   parser.add_argument('-thresh',  type=float, required=True)

   
   args = parser.parse_args()

   get_spec(num=args.num, base=args.base, savedir=args.savedir, res=args.res, numlos=args.numlos, thresh=args.thresh)


