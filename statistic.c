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

#define PBINS 21
int output(double *array, int size, char *suffix, char *outdir);
void help(void);

int main(int argc, char **argv)
{
  double redshift;
  double *tau_H1=NULL;
  double *tau_He2=NULL;
  int  NumLos=0;
  int UsedLos=0;
  FILE *input;
  /*Which statistic to use: 1 is pdf, 2 is power,
   * 3 is transverse power, 4 is bispectrum.
   * Only 1 and 2 are implemented.*/
  int statistic=2;
  int rescale=1;
  double scale=1.0;
/*   double tau_effs[11]={0.178000, 0.2192, 0.2714000, 0.3285330, 0.379867, 0.42900, 0.513000, 0.600400, 0.657800,  0.756733,  0.896000}; */
/*   double tau_eff; */
  double flux_power[(NBINS+1)/2];
  double flux_pdf[PBINS];
  char *inname=NULL;
  char *outdir=NULL;
  char suffix[30];
  char c;
  /*Make sure stdout is line buffered even when not
   * printing to a terminal but, eg, perl*/
  setlinebuf(stdout);
  while((c = getopt(argc, argv, "s:o:i:n:u:rh")) !=-1)
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
        case 'u':
           UsedLos = atoi(optarg);
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
  if(UsedLos == 0)
          UsedLos = NumLos;
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

  if(!(tau_H1 = (double *) calloc((UsedLos * NBINS) , sizeof(double))))
  {
      fprintf(stderr, "Failed to allocate memory!\n");
      exit(1);
  }
#ifdef HELIUM
  if(!(tau_He2 = (double *) calloc((UsedLos * NBINS) , sizeof(double))))
  {
      fprintf(stderr, "Failed to allocate helium memory!\n");
      exit(1);
  }
#endif
  if(!(input=fopen(inname,"rb"))){
        fprintf(stderr, "Could not open file %s for reading!\n",inname);
        exit(2);
  }
  fread(&redshift,sizeof(double),1,input);
  fseek(input,sizeof(double)*NBINS*NumLos*4,SEEK_CUR);
  if(fread(tau_H1,sizeof(double),NBINS*UsedLos,input) != NBINS*UsedLos)     /* HI optical depth */
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
  printf("NumLos=%d tau_H1[0]=%g tau_H1[N]=%g\n",UsedLos,tau_H1[0],tau_H1[NBINS*UsedLos-1]);
  /*Calculate mean flux*/
  /*Changing mean flux by a factor of ten changes the P_F by a factor of three*/
/*   tau_eff=tau_effs[(int)(redshift-2.2)*5]; */
  if(rescale)
  {
    scale=mean_flux(tau_H1, NBINS*UsedLos,exp(-TAU_EFF),1e-5 );
    printf("scale=%g\n",scale);
  }
  /*If no rescale, we output the non-rescaled power spectrum as well*/
  if(statistic == 2){
      calc_power_spectra(flux_power,tau_H1,scale,TAU_EFF,UsedLos);
      sprintf(suffix,"_flux_power.txt");
      if(output(flux_power, (NBINS+1)/2,suffix, outdir))
          exit(1);
  }
  if(statistic == 1){
      calc_pdf(flux_pdf, tau_H1,scale,UsedLos);
      sprintf(suffix,"_flux_pdf.txt");
      if(output(flux_pdf, PBINS,suffix, outdir))
          exit(1);
  }
  /*Average the power spectrum*/
  free(tau_H1);
#ifdef HELIUM
  free(tau_He2);
#endif
  return 0;
}

int output(double *array, int size, char *suffix, char *outdir)
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

void calc_pdf(double *flux_pdf, double *tau_H1, double scale, int NumLos)
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
        flux_pdf[i]*=(PBINS-1);
    }
    return;
}

/*****************************************************************************/
void help()
{
           fprintf(stderr, "Usage: ./statistic -n NUMLOS -i input spectra -o output_directory\n"
                           "-r turns off rescaling -s statistic: 1 is pdf, 2 is power,\n"
                           "3 is transverse power, 4 is bispectrum");
           return;
}


