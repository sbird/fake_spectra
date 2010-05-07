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

#define RAW_SPECTRA 1
#define _GNU_SOURCE
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "global_vars.h"
#include "parameters.h"
#define PTS 30
#define TAU_MEAN_MIN 0.2
#define TAU_MEAN_MAX 0.5
#define TAU_SLOPE_MIN 2.0
#define TAU_SLOPE_MAX 4.0
#define TAU_MIN TAU_MEAN_MIN*pow((1.0+redshift)/4.0,TAU_SLOPE_MIN)
#define TAU_MAX TAU_MEAN_MAX*pow((1.0+redshift)/4.0,TAU_SLOPE_MAX)

int main(int argc, char **argv)
{
  int  NumLos=0;
  FILE *input;
  FILE *output;
  char *inname=NULL;
  char *outdir=NULL;
  char *outname=NULL;
  char c;
  int j,jj,n;
  /*Make sure stdout is line buffered even when not 
   * printing to a terminal but, eg, perl*/
  setlinebuf(stdout);
  while((c = getopt(argc, argv, "o:i:n:h")) !=-1)
  {
    switch(c)
      {
        case 'o':
           outdir=optarg;
           break;
        case 'i':
           inname=optarg;
           break;
        case 'n':
           NumLos=atoi(optarg);
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
  if( !outdir || !inname)
  {
         fprintf(stderr, "Specify output (%s) and input (%s) directories.\n",outdir, inname);
         help();
         exit(99);
  }
  InitLOSMemory(NumLos);
  input=fopen(inname,"rb");
  fread(&redshift,sizeof(double),1,input);
  fseek(input,sizeof(double)*NBINS*NumLos*4,SEEK_CUR);
  fread(tau_H1,sizeof(double),NBINS*NumLos,input);    /* HI optical depth */
#ifdef HELIUM
  fseek(input,sizeof(double)*NBINS*NumLos*3,SEEK_CUR);
  fread(tau_He2,sizeof(double),NBINS*NumLos,input);   /* HeII optical depth */
#endif
  fclose(input);
  printf("NumLos=%d tau_H1[0]=%g tau_H1[N]=%g\n",NumLos,tau_H1[0],tau_H1[NBINS*NumLos-1]);
  /*Calculate mean flux*/
  /*Changing mean flux by a factor of ten changes the P_F by a factor of three*/
  /*Do this even if not rescaling, as we want to know the mean flux*/
  pl=rfftw_create_plan(NBINS,FFTW_REAL_TO_COMPLEX, FFTW_MEASURE | FFTW_THREADSAFE);

  inname=basename(inname);
  outname=malloc((strlen(outdir)+strlen(inname)+80)*sizeof(char));
  for(n=0; n<PTS; n++)
  {
          double scale=1;
          float flux_power_avg[(NBINS+1)/2];
          const double tau_eff=TAU_MIN+(TAU_MAX-TAU_MIN)/PTS*n;
          const double obs_flux=exp(-tau_eff);
          scale=mean_flux(tau_H1, NBINS*NumLos,obs_flux,0.001 );
          printf("scale=%g\n",scale);
          for(j=0; j<(NBINS+1)/2;j++)
            flux_power_avg[j]=0;
          /*If no rescale, we output the non-rescaled power spectrum as well*/
          calc_power_spectra(flux_power,tau_H1,scale,NumLos);
          /*Average the power spectrum*/
          for(j=0; j<(NBINS+1)/2;j++)
          {
              for(jj=0; jj<NumLos; jj++)
                 flux_power_avg[j]+=flux_power[(NBINS+1)/2*jj+j];
              flux_power_avg[j]/=NumLos;
          }
          sprintf(outname,"%s/%s_%g_flux_power.txt",outdir,inname,tau_eff);
          printf("Outputting %s\n",outname);
          output=fopen(outname,"w");
          for(j=0; j<(NBINS+1)/2;j++)
          {
                  /*First value is k*/
             fprintf(output, "%g %g\n", j+0.5, flux_power_avg[j]);
          }
          fclose(output);
  }
  free(outname);
  FreeLOSMemory();
  fftw_destroy_plan(pl);
  return 0;
}

void calc_power_spectra(float *flux_power, double *tau_H1,double scale,int NumLos)
{
    int iproc;
#pragma omp parallel
    {
      /*Perform the scaling*/
      if(scale != 1.0){
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
       return;
}
/*****************************************************************************/
void InitLOSMemory(int NumLos)
{  
  tau_H1       = (double *) calloc((NumLos * NBINS) , sizeof(double));
  flux_power   = (float *) calloc(NumLos * (NBINS+1)/2, sizeof(float));
    if( !tau_H1 || !flux_power  )
  {
      fprintf(stderr, "Failed to allocate memory!\n");
      exit(1);
  }
#ifdef HELIUM 
  tau_He2       = (double *) calloc((NumLos * NBINS) , sizeof(double)); 
  if(! tau_He2 )
  {
      fprintf(stderr, "Failed to allocate helium memory!\n");
      exit(1);
  }
#endif
}
/*****************************************************************************/

/*****************************************************************************/
void FreeLOSMemory(void)
{  
  free(flux_power);
  free(tau_H1   ) ;
#ifdef HELIUM 
  free(tau_He2   );
#endif
}
/*****************************************************************************/
void help()
{
           fprintf(stderr, "Usage: ./extract -n NUMLOS -i filename -o output_directory\n");
           return;
}
