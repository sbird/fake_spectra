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
#include <string.h>
#include <unistd.h>
#include "global_vars.h"
#include "parameters.h"


/* Here we load a snapshot file. It can be distributed
 * onto several files (for files>1) */
/**********************************************************************/
int main(int argc, char **argv)
{
  int Npart, NumLos=0, files=0;
  FILE *output;
  int rescale=1;
  char *outname=NULL;
  char *outdir=NULL;
  char *indir=NULL;
  char c;
#ifndef RAW_SPECTRA
  double obs_flux,scale;
  float flux_power_avg[(NBINS+1)/2];
  int iproc,j,jj;
#endif
  /*Make sure stdout is line buffered even when not 
   * printing to a terminal but, eg, perl*/
  setlinebuf(stdout);
  while((c = getopt(argc, argv, "f:o:i:n:r")) !=-1)
  {
    switch(c)
      {
        case 'f':
           files=atoi(optarg);
           break;
        case 'o':
           outdir=optarg;
           break;
        case 'i':
           indir=optarg;
           break;
        case 'n':
           NumLos=atoi(optarg);
           break;
        case 'r':
           rescale=0;
           break;
        case 'h':
        case '?':
           help();
        default:
           exit(1);
      }
  }
  if(NumLos <=0 || files <=0)
  {
          fprintf(stderr,"Need NUMLOS and NUMFILES >0\n");
          help();
          exit(99);
  
  }
  if( !outdir || !indir)
  {
         fprintf(stderr, "Specify output (%s) and input (%s) directories.\n",outdir, indir);
         help();
         exit(99);
  }
  Npart=load_snapshot(indir, files);
  InitLOSMemory(NumLos);
  if(!PARTTYPE)
    SPH_interpolation(NumLos,Npart);
  /*Free the particle list once we don't need it*/
  free(P);
  /*If the spectra flag is set, output the raw spectra to a file, 
   * instead of the flux power spectrum directly.*/
#ifdef RAW_SPECTRA
  /*Output to a file*/
  outname=malloc((strlen(outdir)+25)*sizeof(char));
  if(!strcpy(outname,outdir) || !(outname=strcat(outname, "_spectra.dat")))
  {
    fprintf(stderr, "Some problem with the strings\n");
    exit(1);
  }
  output=fopen(outname,"w");
  fwrite(&redshift,sizeof(double),1,output);
  fwrite(Delta,sizeof(double),NBINS*NumLos,output);     /* gas overdensity */
  fwrite(n_H1,sizeof(double),NBINS*NumLos,output);      /* n_HI/n_H */
  fwrite(temp_H1,sizeof(double),NBINS*NumLos,output);   /* T [K], HI weighted */
  fwrite(veloc_H1,sizeof(double),NBINS*NumLos,output);  /* v_pec [km s^-1], HI weighted */
  fwrite(tau_H1,sizeof(double),NBINS*NumLos,output);    /* HI optical depth */
#ifdef HELIUM
  fwrite(n_He2,sizeof(double),NBINS*NumLos,output);     /* n_HeII/n_H */
  fwrite(temp_He2,sizeof(double),NBINS*NumLos,output);   /* T [K], HeII weighted */
  fwrite(veloc_He2,sizeof(double),NBINS*NumLos,output); /* v_pec [km s^-1], HeII weighted */
  fwrite(tau_He2,sizeof(double),NBINS*NumLos,output);   /* HeII optical depth */
#endif
  fwrite(posaxis,sizeof(double),NBINS,output);          /* pixel positions, comoving kpc/h */
  fwrite(velaxis,sizeof(double),NBINS,output);          /* pixel positions, km s^-1 */
  fclose(output);
#else /*If outputting flux power directly*/
  /*Calculate mean flux*/
  /*Changing mean flux by a factor of ten changes the P_F by a factor of three*/
  if(rescale){
    obs_flux= exp(-TAU_EFF);
    scale=mean_flux(tau_H1, NBINS*NumLos,obs_flux,0.001 );
    printf("scale=%g\n",scale);
  }
  for(j=0; j<(NBINS+1)/2;j++)
    flux_power_avg[j]=0;
  pl=rfftw_create_plan(NBINS,FFTW_REAL_TO_COMPLEX, FFTW_MEASURE | FFTW_THREADSAFE);
#pragma omp parallel
  {
    if(rescale){
      /*Perform the scaling*/
      #pragma omp for schedule(static, THREAD_ALLOC)
      for(iproc=0; iproc<NBINS*NumLos; iproc++)
      {
        tau_H1[iproc]*=scale;
      }
    }
    /*Calculate power spectrum*/
    #pragma omp for schedule(static, THREAD_ALLOC)
    for(iproc=0;iproc<NumLos;iproc++)
    {
        const double * tau_H1_local = &tau_H1[iproc*NBINS];
        float flux_power_local[(NBINS+1)/2];
        float flux_H1_local[NBINS];
        int i,ii;
        /* Calculate flux and flux power spectrum */
        for(i=0; i<NBINS; i++)
        {
           flux_H1_local[i]=exp(-tau_H1_local[i]);
        }
        powerspectrum(NBINS, flux_H1_local, flux_power_local);
        /*Write powerspectrum*/
        for(i=0; i<(NBINS+1)/2;i++)
        {
            ii=i+(NBINS+1)/2*iproc;
            flux_power[ii]=flux_power_local[i];
        }
    }/*End loop*/
  }/*End parallel block*/
  /*Average the power spectrum*/
  for(j=0; j<(NBINS+1)/2;j++)
  {
      for(jj=0; jj<NumLos; jj++)
         flux_power_avg[j]+=flux_power[(NBINS+1)/2*jj+j];
      flux_power_avg[j]/=NumLos;
  }
  printf("Outputting average flux power spectrum\n");
  outname=malloc((strlen(outdir)+25)*sizeof(char));
  if(!strcpy(outname,outdir) || !(outname=strcat(outname, "_flux_power.txt")))
  {
    fprintf(stderr, "Some problem with the strings\n");
    exit(1);
  }
  output=fopen(outname,"w");
  for(j=0; j<(NBINS+1)/2;j++)
  {
          /*First value is k*/
     fprintf(output, "%g %g\n", j+0.5, flux_power_avg[j]);
  }
  fclose(output);
  fftw_destroy_plan(pl);
  #endif /* RAW_SPECTRA*/
  /*Free other memory*/
  free(outname);
  FreeLOSMemory();
  return 0;
}
/**********************************************************************/

void help()
{
           fprintf(stderr, "Usage: ./extract -f NUMFILES -n NUMLOS -i filename (ie, without the .0) -o output_file (_flux_power.txt or _spectra.dat will be appended)\n"
                           "-r turns off mean flux rescaling\n");
           return;
}
