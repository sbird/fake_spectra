import fake_spectra
#print("fake_spectra imported")
from fake_spectra.randspectra import RandSpectra
#print("haloassigned_spectra imported")
from fake_spectra.plot_spectra import PlottingSpectra
#print("PlottingSpectra imported")
from mpi4py import MPI
import numpy as np

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

size_col_den = 78755 # size of array col_den

### the array to store added column density in
col_den_added = np.empty(size_col_den, 'd')


if rank ==0 :

    for i in range(0,100):
        
        istr = str(i)
        rr = RandSpectra(34, "/rhome/mqezl001/bigdata/TNG/TNG100-1/output", file_num = i, thresh = 0.0, kernel='tophot',ndla = 1000, numlos=1000,savedir="/rhome/mqezl001/bigdata/TNG/TNG100-1/postprocessing/randspectra/Snap_034/", savefile="spectra_34."+istr+".hdf5")
        
        num_not_DLA = 1000 # just for first iteration  of the while below

        while(num_not_DLA > 0) :

            ### Here, get not_DLA_indices and if there are still some, regenerate those indices
            ### Replacing part
            rr.replace_not_DLA(ndla, thresh, elem=elem, ion=ion)
        
            tau_1215 = rr.get_tau("H", 1,1215)
            tau_1025 = rr.get_tau("H", 1, 1025)
            col_den_total = rr.get_col_density("H",1)
            col_den_HI = rr.get_col_density('H', -1)
        
            ### Send col_den to root rank

            comm_Reduce(col_den_HI, col_den_added, op=MPI.SUM, root=2)
        
            ### Recieve not_DLA_indices from manager rank
            
            num_not_DLA = comm.bcast(root=2)
                    
            not_DLA_indices = np.empty(num_not_DLA, dtype='d')
            comm.Bcast(not_DLA_indices, root=2)
            

"""
if rank ==1 :

    for i in range(100,200):

        istr = str(i)
        rr = RandSpectra(34, "/rhome/mqezl001/bigdata/TNG/TNG100-1/output", file_num = i, thresh = 0.0, kernel='tophot',ndla = 1000, numlos=1000,savedir="/rhome/mqezl001/bigdata/TNG/TNG100-1/postprocessing/randspectra/Snap_034/", savefile="spectra_34."+istr+".hdf5")

        
        tau_1215 = rr.get_tau("H",1,1215)
        tau_1025 = rr.get_tau("H", 1, 1025)
        col_den_total = rr.get_col_density("H",1)
        col_den_HI = rr.get_col_density('H', -1)
        
        comm_Send(col_den_HI, dest=2, tag=24)

"""

if rank ==2 :
    

    While (num_not_DLA > 0) :

        col_den_added = np.empty(size_col_den, dtype='d')
        #temp_col_den_added = np.empty(size_col, dtype='d')

        ### Recieve col_den from all other ranks and add them together
        
        comm.Reduce(col_den_added, col_den_added, op = MPI.SUM, root = 2)
        
        ### indices which should be replaced
        not_DLA_indices = np.where(np.sum(col_den_added[:], axis=1) < 10**20.3)
        num_not_DLA = np.size(not_DLA_indices)
        ### Here, Broadcast not_DLA_indices to each rank to be regenrated
        

        comm.bcast(np.size(num_not_DLA, root=2)
        comm.Bcast(not_DLA_indices, root = 2)




    ### Save Final Data

    rr.savefile()
        




