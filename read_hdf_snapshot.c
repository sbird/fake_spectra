#ifdef HDF5

#include <hdf5.h>
#include <hdf5_hl.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "global_vars.h"
#include "parameters.h"

#ifndef N_TYPE
        #define N_TYPE 6
#endif

/*Open a file for reading to check it exists*/
int file_readable(const char * filename)
{   
     FILE * file; 
     if ((file = fopen(filename, "r"))){
          fclose(file);
          return 1;
     }
     return 0;
}

/*Routine that is a wrapper around HDF5's dataset access routines to do error checking. Returns the length on success, 0 on failure.*/
hsize_t get_single_dataset(const char *name, void * data_ptr, hid_t * hdf_group,int fileno){
          int rank;
          hsize_t vlength;
          size_t type_size;
          H5T_class_t class_id;
          if (H5LTget_dataset_ndims(*hdf_group, name, &rank) < 0 || rank != 1){
             fprintf(stderr, "File %d: Rank of %s is %d !=1\n",fileno,name, rank);
             return 0;
          }
          H5LTget_dataset_info(*hdf_group, name, &vlength, &class_id, &type_size);
          if(type_size != 4 || class_id != H5T_FLOAT || H5LTread_dataset_float(*hdf_group, name, data_ptr) < 0 ){
              fprintf(stderr, "File %d: Failed reading %s (%lu)\n",fileno,name, (uint64_t)vlength);
              return 0;
          }
          return vlength;
}

/*A similar wrapper around HDF5's dataset access routines to do error checking. Returns the length on success, 0 on failure.*/
hsize_t get_triple_dataset(const char *name, void * data_ptr, hid_t * hdf_group,int fileno){
          int rank;
          hsize_t vlength[2];
          size_t type_size;
          H5T_class_t class_id;
          if (H5LTget_dataset_ndims(*hdf_group, name, &rank) < 0 || rank != 2){
             fprintf(stderr, "File %d: Rank of %s is %d !=2\n",fileno,name, rank);
             return 0;
          }
          H5LTget_dataset_info(*hdf_group, name, &vlength[0], &class_id, &type_size);
          if(type_size != 4 || class_id != H5T_FLOAT || vlength[1] != 3 || H5LTread_dataset_float(*hdf_group, name, data_ptr) < 0 ){
              fprintf(stderr, "File %d: Failed reading %s (%lu)\n",fileno,name, (uint64_t)vlength[0]);
              return 0;
          }
          return vlength[0];
}


/* this routine loads particle data from an HDF5 snapshot. 
 * (A snapshot may be distributed
 * into multiple files. */
int load_hdf5_snapshot(char *infname, pdata *P,
  double  *atime, double *redshift, double * Hz, double *box100, double *h100, double *omegab)
{
  int i, fileno=0;
  int64_t NumPart,NumRead=0;
  int npart[N_TYPE];
  double mass[N_TYPE];
  char name[16];
  int flag_cooling;
  int64_t npart_all[N_TYPE];
  double Omega0, OmegaLambda;
  hid_t hdf_group,hdf_file;
  hsize_t length;
  /*Copy of input filename for extension*/
  char fname[strlen(infname)+10];
  char ffname[strlen(infname)+16];
  /*Switch off error handling so that we can check whether a
   * file is HDF5 */
  /* Save old error handler */
  hid_t error_stack=0;
  herr_t (*old_func)(hid_t, void*);
  void *old_client_data;
  H5Eget_auto(error_stack, &old_func, &old_client_data);
  /* Turn off error handling */
  H5Eset_auto(error_stack, NULL, NULL);

  /*Were we handed an HDF5 file?*/
  if(H5Fis_hdf5(infname) <= 0){
     /*If we weren't, were we handed an HDF5 file without the suffix?*/
     strncpy(fname, infname,strlen(infname));
     strncpy(fname+strlen(infname), ".0.hdf5\0",10);
     if (H5Fis_hdf5(fname) <= 0)
        return -1;
  }
  else{
     strncpy(fname, infname,strlen(infname));
  }

  /* Restore previous error handler */
  H5Eset_auto(error_stack, old_func, old_client_data);

  /*See if we have been handed the first file of a set: 
   * our method for dealing with this closely mirrors 
   * HDF5s family mode, but we cannot use this, because
   * our files may not all be the same size.*/
  char *zero = strstr(fname,".0.hdf5");
  /*Replace a possible 0.hdf5 in the filename 
   * with a printf style specifier for opening*/
  if(zero)
    strncpy(zero, ".%d.hdf5\0",strlen(zero)+3);

  /*First open first file to get header properties*/ 
  sprintf(ffname,fname,fileno);
  hdf_file=H5Fopen(ffname,H5F_ACC_RDONLY,H5P_DEFAULT);
  if(hdf_file < 0){
        return -1;
  }
  if ( (hdf_group=H5Gopen(hdf_file,"/Header",H5P_DEFAULT)) < 0) {
        H5Fclose(hdf_file);
        return -1;
  }
  /* Read some header functions */
  
  if(H5LTget_attribute_double(hdf_group,".","Time",atime) ||
     H5LTget_attribute_double(hdf_group,".","Redshift", redshift) ||
     H5LTget_attribute_double(hdf_group,".","BoxSize", box100) ||
     H5LTget_attribute_double(hdf_group,".","HubbleParam", h100) ||
     H5LTget_attribute_double(hdf_group,".","Omega0", &Omega0) ||
     H5LTget_attribute_double(hdf_group,".","OmegaLambda", &OmegaLambda) ||
     H5LTget_attribute_int(hdf_group,".","Flag_Cooling",&flag_cooling)){
          fprintf(stderr,"Failed to read some header value\n");
      H5Gclose(hdf_group);
      H5Fclose(hdf_file);
      return -1;
  }
  (*Hz)=100.0*(*h100) * sqrt(1.+Omega0*(1./(*atime)-1.)+OmegaLambda*((pow(*atime,2)) -1.))/(*atime);
  /*Get the total number of particles*/
  H5LTget_attribute(hdf_group,".","NumPart_Total",H5T_NATIVE_INT, &npart);
  for(i = 0; i< N_TYPE; i++)
          npart_all[i]=npart[i];
  H5LTget_attribute(hdf_group,".","NumPart_Total_HighWord",H5T_NATIVE_INT, &npart);
  for(i = 0; i< N_TYPE; i++)
          npart_all[i]+=(1L<<32)*npart[i];
  H5LTget_attribute(hdf_group,".","MassTable",H5T_NATIVE_DOUBLE, mass);
  
  /*Close header*/
  H5Gclose(hdf_group);
  H5Fclose(hdf_file);
  
  NumPart = npart_all[PARTTYPE];
  if(NumPart <= 0)
          return 0;
  printf("NumPart=[%ld,%ld,%ld,%ld,%ld,%ld], ",npart_all[0],npart_all[1],npart_all[2],npart_all[3],npart_all[4],npart_all[5]);
  printf("Masses=[%g %g %g %g %g %g], ",mass[0],mass[1],mass[2],mass[3],mass[4],mass[5]);
  printf("Redshift=%g, Ω_M=%g Ω_L=%g\n",(*redshift),Omega0,OmegaLambda);
  printf("Expansion factor = %f\n",(*atime));
  printf("Hubble = %g Box=%g \n",(*h100),(*box100));
  
  if(!(alloc_parts(P,NumPart)))
  {
    fprintf(stderr,"Failed to allocate memory.\n\n");
    exit(1);
  }
  /*Loop over files. Keep going until we run out, skipping over broken files.
   * The call to file_readable is an easy way to shut up HDF5's error message.*/
  while(file_readable(ffname) && H5Fis_hdf5(ffname) > 0){
      hdf_file=H5Fopen(ffname,H5F_ACC_RDONLY,H5P_DEFAULT);
      if(hdf_file < 0){
            fileno++;
            sprintf(ffname,fname,fileno);
            continue;
      }
      /*Open particle data*/
      snprintf(name,16,"/PartType%d",PARTTYPE);

      if ( (hdf_group=H5Gopen(hdf_file,name,H5P_DEFAULT)) < 0) {
            fileno++;
            sprintf(ffname,fname,fileno);
            H5Fclose(hdf_file);
            continue;
      }
      
      /* Read position and velocity*/
      length = get_triple_dataset("Coordinates",(*P).Pos+3*NumRead,&hdf_group,fileno);
      if(length == 0)
              goto exit;
      printf("Reading File %d (%lu particles)\n", fileno,(uint64_t)length);
      if(length != get_triple_dataset("Velocities",(*P).Vel+3*NumRead,&hdf_group,fileno))
              goto exit;
      
      /* Particle masses  */
      if(mass[PARTTYPE])
            for(i=0; i< length;i++)
               (*P).Mass[i] = mass[PARTTYPE];
      else
         if (length != get_single_dataset("Masses",(*P).Mass+NumRead,&hdf_group,fileno))
                 goto exit;
      (*omegab) = (*P).Mass[0]/((*P).Mass[0]+mass[1])*Omega0;
      /*Seek past the last masses*/
      if(PARTTYPE == 0)
        { 
          /*The internal energy of all the Sph particles is read in */
         if (length != get_single_dataset("InternalEnergy",(*P).U+NumRead,&hdf_group,fileno))
                 goto exit;
          /* The free electron fraction */
          if(flag_cooling)
            {
              /* Some versions of Gadget have Ne, some have NHP, NHEP and NHEPP*/
              /* If this were supported, I would use that the universe is neutral, so 
               * NE = NHP + NHEP +2 NHEPP*/
              if (length != get_single_dataset("ElectronAbundance",(*P).Ne+NumRead,&hdf_group,fileno))
                 goto exit;
              /* The HI fraction, nHI/nH */
              if (length != get_single_dataset("NeutralHydrogenAbundance",(*P).NH0+NumRead,&hdf_group,fileno))
                 goto exit;
            }
         /* The smoothing length */
         if (length != get_single_dataset("SmoothingLength",(*P).h+NumRead,&hdf_group,fileno))
            goto exit;
        }
      NumRead += length;
exit:
      H5Gclose(hdf_group);
      H5Fclose(hdf_file);
      fileno++;
      sprintf(ffname,fname,fileno);
  }
  printf("P[%d].Pos = [%g %g %g]\n", 0, (*P).Pos[0], (*P).Pos[1],(*P).Pos[2]);
  printf("P[%d].Vel = [%g %g %g]\n", 0, (*P).Vel[0], (*P).Vel[1],(*P).Vel[2]);
  printf("P[%ld].Mass = %e Ω_B=%g\n\n", NumRead, (*P).Mass[0],(*omegab));
  printf("P[%ld].U = %f\n\n", NumRead, (*P).U[NumRead-1]);
  printf("P[%ld].Ne = %e\n", NumRead, (*P).Ne[NumRead-1]);
  printf("P[%ld].NH0 = %e\n", NumRead, (*P).NH0[NumRead-1]);
  printf("P[%ld].h = %f\n",NumRead, (*P).h[NumRead-1]);
  return NumRead;
}
#endif //HDF5
