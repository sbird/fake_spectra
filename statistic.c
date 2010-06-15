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

#define PBINS 20
int output(float *array, int size, char *suffix, char *outdir);

int main(int argc, char **argv)
{
  int  NumLos=0;
  FILE *input;
  /*Which statistic to use: 1 is pdf, 2 is power,
   * 3 is transverse power, 4 is bispectrum.
   * Only 1 and 2 are implemented.*/
  int statistic=2;
  int rescale=1;
  double scale=1.0;
  float flux_power[(NBINS+1)/2];
  float flux_pdf[PBINS];
  char *inname=NULL;
  char *outdir=NULL;
  char suffix[30];
  char c;
  /*Make sure stdout is line buffered even when not
   * printing to a terminal but, eg, perl*/
  setlinebuf(stdout);
  while((c = getopt(argc, argv, "s:o:i:n:rh")) !=-1)
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
        case 's':
           statistic=atoi(optarg);
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
  if(NumLos <=0){
          fprintf(stderr,"Need NUMLOS >0\n");
          help();
          exit(99);
  }
  if( !outdir || !inname){
         fprintf(stderr, "Specify output (%s) and input (%s) directories.\n",outdir, inname);
         help();
         exit(99);
  }
  if(statistic !=1 && statistic !=2){
          fprintf(stderr, "Only flux pdf and flux power are implemented.\n");
          help();
          exit(99);
  }

  InitLOSMemory(NumLos);
  if(!(input=fopen(inname,"rb"))){
        fprintf(stderr, "Could not open file %s for reading!\n",inname);
        exit(2);
  }
  fread(&redshift,sizeof(double),1,input);
  fseek(input,sizeof(double)*NBINS*NumLos*4,SEEK_CUR);
  if(fread(tau_H1,sizeof(double),NBINS*NumLos,input) != NBINS*NumLos)     /* HI optical depth */
  {
          fprintf(stderr, "Could not read spectra!\n");
          exit(2);
  }
#ifdef HELIUM
  fseek(input,sizeof(double)*NBINS*NumLos*3,SEEK_CUR);
  if(fread(tau_He2,sizeof(double),NBINS*NumLos,input) !=NBINS*NumLos)   /* HeII optical depth */
  {
          fprintf(stderr, "Could not read spectra!\n");
          exit(2);
  }
#endif
  fclose(input);
  printf("NumLos=%d tau_H1[0]=%g tau_H1[N]=%g\n",NumLos,tau_H1[0],tau_H1[NBINS*NumLos-1]);
  /*Calculate mean flux*/
  /*Changing mean flux by a factor of ten changes the P_F by a factor of three*/
  if(rescale)
  {
    scale=mean_flux(tau_H1, NBINS*NumLos,exp(-TAU_EFF),0.001 );
    printf("scale=%g\n",scale);
  }
  /*If no rescale, we output the non-rescaled power spectrum as well*/
  if(statistic == 2){
      pl=rfftw_create_plan(NBINS,FFTW_REAL_TO_COMPLEX, FFTW_MEASURE | FFTW_THREADSAFE);
      calc_power_spectra(flux_power,tau_H1,scale,NumLos);
      fftw_destroy_plan(pl);
      if(rescale)
              sprintf(suffix,"_flux_power.txt");
      else
              sprintf(suffix,"_no_rescale_flux_power.txt");
      if(output(flux_power, (NBINS+1)/2,suffix, outdir))
          exit(1);
  }
  if(statistic == 1){
      calc_pdf(flux_pdf, tau_H1,scale,NumLos);
      if(rescale)
              sprintf(suffix,"_flux_pdf.txt");
      else
              sprintf(suffix,"_no_rescale_flux_pdf.txt");
      if(output(flux_pdf, PBINS,suffix, outdir))
          exit(1);
  }
  /*Average the power spectrum*/
  FreeLOSMemory();
  return 0;
}

int output(float *array, int size, char *suffix, char *outdir)
{
     char *outname=NULL;
     FILE *out;
     int j;
     outname=malloc((strlen(outdir)+strlen(suffix)+2)*sizeof(char));
     if(!strcpy(outname,outdir) || !(outname=strcat(outname, suffix)))
     {
       fprintf(stderr, "Some problem with the strings\n");
       return 1;
     }
     if(!(out=fopen(outname,"w")))
     {
       fprintf(stderr, "Error opening %s\n",outname);
       return 1;
     }
     /*First value is k*/
     for(j=0; j<size;j++)
        fprintf(out, "%g %g\n", j+0.5, array[j]);
     fclose(out);
     return 0;
}

void calc_power_spectra(float *flux_power, double *tau_H1,double scale,int NumLos)
{
    int iproc;
    for(iproc=0; iproc<(NBINS+1)/2;iproc++)
      flux_power[iproc]=0;
#pragma omp parallel
    {
      /*Perform the scaling*/
       #pragma omp for schedule(static, THREAD_ALLOC)
       for(iproc=0; iproc<NBINS*NumLos; iproc++)
       {
         tau_H1[iproc]*=scale;
       }
       /*Calculate power spectrum*/
       #pragma omp for schedule(static, THREAD_ALLOC)
       for(iproc=0;iproc<NumLos;iproc++)
       {
           const double * tau_H1_local = &tau_H1[iproc*NBINS];
           float flux_power_local[(NBINS+1)/2];
           float flux_H1_local[NBINS];
           int i;
           /* Calculate flux and flux power spectrum */
           for(i=0; i<NBINS; i++)
           {
              flux_H1_local[i]=exp(-tau_H1_local[i]);
           }
           powerspectrum(NBINS, flux_H1_local, flux_power_local);
           /*Write powerspectrum*/
           for(i=0; i<(NBINS+1)/2;i++)
           {
               #pragma omp atomic
               flux_power[i]+=flux_power_local[i];
           }
       }/*End loop*/
     }/*End parallel block*/
     /*Average the power spectrum*/
     for(iproc=0; iproc<(NBINS+1)/2;iproc++)
     {
         flux_power[iproc]/=NumLos;
     }
     return;
}

void calc_pdf(float *flux_pdf, double *tau_H1, double scale, int NumLos)
{
    int i;
    for(i=0;i<PBINS; i++)
        flux_pdf[i]=0;
    /* Calculate flux pdf */
    for(i=0;i<NBINS*NumLos;i++)
        flux_pdf[(int)round(exp(-scale*tau_H1[i])*(PBINS-1))]++;
    /*Normalise*/
    for(i=0;i<PBINS;i++){
        flux_pdf[i]/=(NumLos*NBINS);
        flux_pdf[i]*=PBINS;
    }
    return;
}
/*****************************************************************************/
void InitLOSMemory(int NumLos)
{
  if(!(tau_H1 = (double *) calloc((NumLos * NBINS) , sizeof(double))))
  {
      fprintf(stderr, "Failed to allocate memory!\n");
      exit(1);
  }
#ifdef HELIUM
  if(!(tau_He2 = (double *) calloc((NumLos * NBINS) , sizeof(double))))
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
  free(tau_H1);
#ifdef HELIUM
  free(tau_He2);
#endif
}
/*****************************************************************************/
void help()
{
           fprintf(stderr, "Usage: ./statistic -n NUMLOS -i input spectra -o output_directory\n"
                           "-r turns off rescaling -s statistic: 1 is pdf, 2 is power,\n"
                           "3 is transverse power, 4 is bispectrum");
           return;
}


