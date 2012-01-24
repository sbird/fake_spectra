/* Copyright (c) 2009, Simeon Bird <spb41@cam.ac.uk>
 *               Based on code (c) 2005 by J. Bolton
 *
 * Permission to use, copy, modify, and/or distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
 * ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 * ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 * OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE. */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <unistd.h>
#include "global_vars.h"
#include "parameters.h"
#ifdef HDF5
  #include <hdf5.h>
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


int main(int argc, char **argv)
{
  int64_t Npart;
  int NumLos=0;

  FILE *output;
  los *los_table=NULL;
  char *ext_table=NULL;
  char *outname=NULL;
  char *outdir=NULL;
  char *indir=NULL;
#ifdef HDF5
  char *fname=NULL;
#endif
  char c;
  int i;
#ifndef NO_HEADER
  int pad[32]={0};
#endif
  double  atime, redshift, Hz, box100, h100, omegab;
  struct particle_data P;
  double * rhoker_H=NULL;
  double * tau_H1=NULL;
  interp H1;
#ifdef HELIUM
  double *tau_He2=NULL;
  interp He2;
#endif
  /*Make sure stdout is line buffered even when not 
   * printing to a terminal but, eg, perl*/
  setlinebuf(stdout);
  while((c = getopt(argc, argv, "o:i:t:n:h")) !=-1)
  {
    switch(c)
      {
        case 'o':
           outdir=optarg;
           break;
        case 'i':
           indir=optarg;
           break;
        case 'n':
           NumLos=atoi(optarg);
           break;
        case 't':
           ext_table=optarg;
           break;
        case 'h':
        case '?':
           help();
        default:
           exit(1);
      }
  }
  if(NumLos <=0)
  {
          fprintf(stderr,"Need NUMLOS >0\n");
          help();
          exit(99);
  
  }
  if( !outdir || !indir)
  {
         fprintf(stderr, "Specify output (%s) and input (%s) directories.\n",outdir, indir);
         help();
         exit(99);
  }
  los_table=malloc(NumLos*sizeof(los));
  if(!los_table){
          fprintf(stderr, "Error allocating memory for sightline table\n");
          exit(2);
  }
  if(InitLOSMemory(&H1, NumLos) || 
      !(rhoker_H = (double *) calloc((NumLos * NBINS) , sizeof(double)))){
                  fprintf(stderr, "Error allocating LOS memory\n");
                  exit(2);
  }
#ifdef HELIUM
  if(InitLOSMemory(&He2, NumLos)){
                  fprintf(stderr, "Error allocating LOS memory\n");
                  exit(2);
  }
#endif
  #ifdef HDF5
    if(!(fname= malloc((strlen(indir)+10)*sizeof(char)))){
        fprintf(stderr, "Failed to allocate string mem\n");
        exit(1);
    }
    /*First open first file to get header properties*/
    if ( find_first_hdf_file(indir,fname) == 0
        && load_hdf5_header(fname, &atime, &redshift, &Hz, &box100, &h100) == 0 ){
        int fileno=0;
        /*Copy of input filename for extension*/
        char ffname[strlen(indir)+16];
        /*See if we have been handed the first file of a set:
         * our method for dealing with this closely mirrors
         * HDF5s family mode, but we cannot use this, because
         * our files may not all be the same size.*/
        char *zero = strstr(fname,".0.hdf5");
        /*Replace a possible 0.hdf5 in the filename
         * with a printf style specifier for opening*/
        if(zero)
          strncpy(zero, ".%d.hdf5\0",strlen(zero)+3);

        populate_los_table(los_table,NumLos, ext_table, box100);
        /*Loop over files. Keep going until we run out, skipping over broken files.
         * The call to file_readable is an easy way to shut up HDF5's error message.*/
        sprintf(ffname,fname,fileno);
        while(file_readable(ffname) && H5Fis_hdf5(ffname) > 0){
              /* P is allocated inside load_hdf5_snapshot*/
              Npart=load_hdf5_snapshot(ffname, &P,&omegab,fileno);
              if(Npart > 0){
           #ifndef HELIUM
              SPH_Interpolation(rhoker_H,&H1,Npart, NumLos,box100, los_table, &P);
           #else
              SPH_Interpolation(rhoker_H,&H1, &He2, Npart, NumLos,box100, los_table, &P);
           #endif
              }
              /*Free the particle list once we don't need it*/
              if(Npart >= 0)
              free_parts(&P);
              fileno++;
              sprintf(ffname,fname,fileno);
        }
    }
    else{
  #endif
        do{
        Npart=load_snapshot(indir, 0,0,&P,&atime, &redshift, &Hz, &box100, &h100, &omegab);
        if(Npart <=0){
                fprintf(stderr,"No data loaded\n");
                exit(99);
        }
        populate_los_table(los_table,NumLos, ext_table, box100);
        /*Do the hard SPH interpolation*/
        #ifndef HELIUM
          SPH_Interpolation(rhoker_H,&H1,Npart, NumLos,box100, los_table, &P);
        #else
          SPH_Interpolation(rhoker_H,&H1, &He2, Npart, NumLos,box100, los_table, &P);
        #endif
        /*Free the particle list once we don't need it*/
        free_parts(&P);
        while
#ifdef HDF5
    }
#endif
  free(los_table);
  if(!(tau_H1 = (double *) calloc((NumLos * NBINS) , sizeof(double)))
  #ifdef HELIUM
   || !(tau_He2 = (double *) calloc((NumLos * NBINS) , sizeof(double)))
  #endif
                  ){
                  fprintf(stderr, "Error allocating memory for tau\n");
                  exit(2);
  }
  printf("Done interpolating, now calculating absorption\n");
#pragma omp parallel
  {
     #pragma omp for
     for(i=0; i<NumLos; i++){
       /*Make a bunch of pointers to the quantities for THIS LOS*/
       interp H1_i=H1;
       #ifdef HELIUM
          interp He2_i=He2;
       #endif
       H1_i.rho+=(i*NBINS);
       H1_i.temp+=(i*NBINS);
       H1_i.veloc+=(i*NBINS);
       #ifndef HELIUM
         Compute_Absorption(tau_H1+(i*NBINS), rhoker_H+(i*NBINS), &H1_i, Hz,h100, box100,atime,omegab);
       #else
         He2_i.rho+=(i*NBINS);
         He2_i.temp+=(i*NBINS);
         He2_i.veloc+=(i*NBINS);
         Compute_Absorption(tau_H1+(i*NBINS), rhoker_H+(i*NBINS), &H1_i, tau_He2+(i*NBINS),&He2_i,Hz,h100,box100,atime,omegab);
       #endif
     }
  }
  /*Output the raw spectra to a file*/ 
  if(!(outname=malloc((strlen(outdir)+29)*sizeof(char))) || !strcpy(outname,outdir) || !(outname=strcat(outname, "_spectra.dat")))
  {
    fprintf(stderr, "Some problem with file output strings\n");
    exit(1);
  }
  if(!(output=fopen(outname,"w")))
  {
          fprintf(stderr, "Error opening %s: %s\n",outname, strerror(errno));
  }
#ifndef NO_HEADER
  /*Write a bit of a header. */
  i=NBINS;
  fwrite(&redshift,sizeof(double),1,output);
  fwrite(&box100,sizeof(double),1,output);
  fwrite(&i,sizeof(int),1,output);
  fwrite(&NumLos,sizeof(int),1,output);
  /*Write some space for future header data: total header size is
   * 128 bytes, with 24 full.*/
  fwrite(&pad,sizeof(int),32-6,output);
#endif
  fwrite(rhoker_H,sizeof(double),NBINS*NumLos,output);     /* gas overdensity */
  if(WriteLOSData(&H1, tau_H1,NumLos, output)
#ifdef HELIUM
    || WriteLOSData(&He2,tau_He2,NumLos, output)
#endif
    ){ 
     fprintf(stderr, "Error writing spectra to disk!\n");
  }
  fclose(output);
  /*Free other memory*/
  free(outname);
  free(tau_H1);
  free(rhoker_H);
  FreeLOSMemory(&H1);
#ifdef HELIUM
  free(tau_He2);
  FreeLOSMemory(&He2);
#endif
  return 0;
}
/**********************************************************************/

void help()
{
           fprintf(stderr, "Usage: ./extract -n NUMLOS -i filename (ie, without the .0) -o output_file (_flux_power.txt or _spectra.dat will be appended)\n"
                  "-t table_file will read line of sight coordinates from a table.\n");
           return;
}

/* Populate the line of sight table, either by random numbers or with some external input. */
void populate_los_table(los * los_table, int NumLos, char * ext_table, double box)
{
        FILE * fh;
        int lines=0;
        int axis;
        float xx, yy, zz;
        /*If we have a file path, load the sightline table from there*/
        if(ext_table){
                const double boxm=box/1000.0;
                if(!(fh=fopen(ext_table, "r")))
                {
                        fprintf(stderr, "Error opening %s: %s\n",ext_table, strerror(errno));
                        exit(3);
                }
                while(lines < NumLos){
                        if(EOF == fscanf(fh, "%d %f %f %f\n", &axis, &xx, &yy, &zz)){
                                fprintf(stderr, "Error reading table: %s. Possibly file is truncated?\n",strerror(errno));
                                exit(3);
                        }
                        if(axis > 3 || axis <0){
                                fprintf(stderr, "Line %d of gives axis %d, which is silly.\n", lines+1, axis);
                                exit(3);
                        }
                        if (xx > boxm || xx < 0 ||
                           yy > boxm || yy < 0 || zz > boxm || zz <0 ){
                                fprintf(stderr, "Line %d of LOS table is: %d %f %f %f, which is silly for boxm %f.\n", lines+1, axis, xx, yy, zz, boxm);
                                exit(3);
                        }
                        los_table[lines].axis=axis;
                        los_table[lines].xx=xx*1000;
                        los_table[lines].yy=yy*1000;
                        los_table[lines].zz=zz*1000;
                        lines++;
                }
        }
        else{
                srand48(241008); /* random seed generator */
                for(lines=0; lines<NumLos; lines++)
                {
                        do	
                        	axis = (int)(drand48()*4);
                        while (axis == 0 || axis==4); 
                        los_table[lines].axis=axis;
                        los_table[lines].xx=drand48()*box;
                        los_table[lines].yy=drand48()*box;
                        los_table[lines].zz=drand48()*box;
                }
        }
        return;
}

/*****************************************************************************/
int InitLOSMemory(interp* species,int NumLos)
{  
  (*species).rho        = (double *) calloc((NumLos * NBINS) , sizeof(double));
  (*species).veloc        = (double *) calloc((NumLos * NBINS) , sizeof(double));
  (*species).temp   = (double *) calloc((NumLos * NBINS) , sizeof(double));
  if(!(*species).rho || !(*species).veloc || !(*species).temp)
      return 1;
  return 0;
}
/*****************************************************************************/

int WriteLOSData(interp* species,double * tau, int NumLos,FILE * output)
{
  int items=0;
  items+=fwrite((*species).rho,sizeof(double),NBINS*NumLos,output);      /* n_HI/n_H */
  items+=fwrite((*species).temp,sizeof(double),NBINS*NumLos,output);   /* T [K], HI weighted */
  items+=fwrite((*species).veloc,sizeof(double),NBINS*NumLos,output);  /* v_pec [km s^-1], HI weighted */
  items+=fwrite(tau,sizeof(double),NBINS*NumLos,output);    /* HI optical depth */
  if(items !=4*NBINS*NumLos)
          return 1;
  return 0;
}

/*****************************************************************************/
void FreeLOSMemory(interp * species)
{  
  free((*species).rho);
  free((*species).veloc);
  free((*species).temp);
}
