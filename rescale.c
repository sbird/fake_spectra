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

#define _GNU_SOURCE
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "statistic.h"
#include "parameters.h"
#define PTS 30
#define TAU_MEAN_MIN 0.2
#define TAU_MEAN_MAX 0.5
#define TAU_SLOPE_MIN 2.0
#define TAU_SLOPE_MAX 4.0
#define TAU_MIN TAU_MEAN_MIN*pow((1.0+redshift)/4.0,TAU_SLOPE_MIN)
#define TAU_MAX TAU_MEAN_MAX*pow((1.0+redshift)/4.0,TAU_SLOPE_MAX)
void help(void);

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
  double redshift;
  double * tau_H1=NULL;
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
  tau_H1       = (double *) calloc((NumLos * NBINS) , sizeof(double));
    if( !tau_H1)
  {
      fprintf(stderr, "Failed to allocate memory!\n");
      exit(1);
  }
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
  inname=basename(inname);
  outname=malloc((strlen(outdir)+strlen(inname)+80)*sizeof(char));
  for(n=0; n<PTS; n++)
  {
          double scale=1.;
          double flux_power_avg[(NBINS+1)/2];
          const double tau_eff=TAU_MIN+n*(TAU_MAX-TAU_MIN)/(double)PTS;
          const double obs_flux=exp(-tau_eff);
          /*Calculate mean flux*/
          scale=mean_flux(tau_H1, NBINS*NumLos,obs_flux,1e-5 );
          printf("tau_eff = %g, scale=%g\n",tau_eff,scale);
          /*If no rescale, we output the non-rescaled power spectrum as well*/
          calc_power_spectra(flux_power_avg,tau_H1,scale,tau_eff,NumLos);
          sprintf(outname,"%s/%d/%s_flux_power.txt",outdir,n, inname);
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
  free(tau_H1   ) ;
  return 0;
}
void help()
{
           fprintf(stderr, "Usage: ./extract -n NUMLOS -i filename -o output_directory\n");
           return;
}
