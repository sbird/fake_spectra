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

/*The value from 0711.1862 is (0.0023±0.0007) (1+z)^(3.65±0.21)*/
#define TAU_EFF 0.0023*pow(1.0+redshift,3.65)

#define PBINS 20
int output(double *array, int size, char *suffix, char *outdir);
void help(void);

int main(int argc, char **argv)
{
  double redshift;
  double *tau_H1=NULL;
#ifdef HELIUM
  double *tau_He2=NULL;
#endif
#ifndef NO_HEADER
  int pad[32];
#endif
  int nbins=NBINS;
  int  NumLos=0;
  int UsedLos=0;
  double box100;
  FILE *input;
  /*Which statistic to use: 1 is pdf, 2 is power,
   * 3 is transverse power, 4 is bispectrum.
   * Only 1 and 2 are implemented.*/
  int statistic=0;
  int rescale=1;
  double scale=1.0;
/*   double tau_effs[11]={0.178000, 0.2192, 0.2714000, 0.3285330, 0.379867, 0.42900, 0.513000, 0.600400, 0.657800,  0.756733,  0.896000}; */
/*   double tau_eff; */
  double *flux_power;
  double flux_pdf[PBINS];
  char *inname=NULL;
  char *outdir=NULL;
  char suffix[30];
  char c;
  /*Make sure stdout is line buffered even when not
   * printing to a terminal but, eg, perl*/
  setlinebuf(stdout);
#ifdef NO_HEADER
  while((c = getopt(argc, argv, "s:o:i:n:u:rh")) !=-1)
#else
  while((c = getopt(argc, argv, "s:o:i:u:rh")) !=-1)
#endif
  {
    switch(c)
      {
        case 'o':
           outdir=optarg;
           break;
        case 'i':
           inname=optarg;
           break;
#ifdef NO_HEADER
        case 'n':
           NumLos=atoi(optarg);
           break;
#endif
        case 's':
           statistic+=atoi(optarg);
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
  if(!inname){
         fprintf(stderr, "Specify input (%s) directory with -i.\n",inname);
         help();
         exit(99);
  }
  if(!outdir){
      outdir = inname;
  }
  if(!(statistic & 1 || statistic & 2)){
          fprintf(stderr, "Only flux pdf and flux power are implemented.\n");
          help();
          exit(99);
  }
  if(!(input=fopen(inname,"rb"))){
        fprintf(stderr, "Could not open file %s for reading!\n",inname);
        exit(2);
  }
  /*read the header*/
  fread(&redshift,sizeof(double),1,input);
#ifndef NO_HEADER
  /*Read a bit of a header. */
  fread(&box100,sizeof(double),1,input);
  fread(&nbins,sizeof(int),1,input);
  fread(&NumLos,sizeof(int),1,input);
  /*Write some space for future header data: total header size is
   * 128 bytes, with 24 full.*/
  fread(&pad,sizeof(int),32-6,input);
#endif
  if(UsedLos == 0)
          UsedLos = NumLos;
  printf("%s using %d sightlines\n",inname,UsedLos);
  if(NumLos <=0){
          fprintf(stderr,"Need NUMLOS >0\n");
          help();
          exit(99);
  }
  if(!(tau_H1 = (double *) calloc((UsedLos * nbins) , sizeof(double)))){
      fprintf(stderr, "Failed to allocate memory!\n");
      exit(1);
  }
  if(!(flux_power = (double *) calloc(((nbins+1)/2) , sizeof(double)))){
      fprintf(stderr, "Failed to allocate memory for flux_power\n");
      exit(1);
  }
#ifdef HELIUM
  if(!(tau_He2 = (double *) calloc((UsedLos * nbins) , sizeof(double))))
  {
      fprintf(stderr, "Failed to allocate helium memory!\n");
      exit(1);
  }
#endif
/*   fseek(input,sizeof(double)*nbins*NumLos*2,SEEK_CUR); */
  if(fread(tau_H1,sizeof(double),nbins*UsedLos,input) != nbins*UsedLos)     /* HI optical depth */
  {
          fprintf(stderr, "Could not read spectra!\n");
          exit(2);
  }
#ifdef HELIUM
  fseek(input,sizeof(double)*nbins*NumLos*3,SEEK_CUR);
  if(fread(tau_He2,sizeof(double),nbins*NumLos,input) !=nbins*NumLos)   /* HeII optical depth */
  {
          fprintf(stderr, "Could not read spectra!\n");
          exit(2);
  }
#endif
  fclose(input);
  printf("NumLos=%d tau_H1[0]=%g tau_H1[N]=%g\n",UsedLos,tau_H1[0],tau_H1[nbins*UsedLos-1]);
  /*Calculate mean flux*/
  /*Changing mean flux by a factor of ten changes the P_F by a factor of three*/
/*   tau_eff=tau_effs[(int)(redshift-2.2)*5]; */

  if(rescale)
  {
    scale=mean_flux(tau_H1, nbins*UsedLos,exp(-TAU_EFF),1e-5 );
    printf("scale=%g\n",scale);
  }
  else
  {
    double mean_flux = 0;
    int i;
    for(i=0; i< nbins*UsedLos; i++)
    {
        mean_flux+=exp(-tau_H1[i]);
    }
    mean_flux/=nbins*UsedLos;
    printf("Mean flux is %g\n",mean_flux);
  }
  /*If no rescale, we output the non-rescaled power spectrum as well*/
  if(statistic & 2){
      calc_power_spectra(flux_power,tau_H1,scale,TAU_EFF,UsedLos);
      sprintf(suffix,"_flux_power.txt");
      if(output(flux_power, (nbins+1)/2,suffix, outdir))
          exit(1);
  }
  if(statistic & 1){
      calc_pdf(flux_pdf, tau_H1,scale,UsedLos,nbins);
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
     char *outname=NULL, *replaceable;
     FILE *out;
     int j;
     outname=malloc((strlen(outdir)+strlen(suffix)+2)*sizeof(char));
     if(!strcpy(outname,outdir) || !(outname=strcat(outname, suffix)))
     {
       fprintf(stderr, "Some problem with the strings\n");
       return 1;
     }
     /*Replace _spectra.dat in the filename with the suffix, instead of appending*/
     replaceable = strstr(outname, "_spectra.dat");
     if(replaceable){
       strcpy(replaceable,suffix);
       *(replaceable+strlen(suffix))='\0';
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

void calc_pdf(double *flux_pdf, double *tau_H1, double scale, int NumLos, int nbins)
{
    int i;
    for(i=0;i<PBINS; i++)
        flux_pdf[i]=0;
    /* Calculate flux pdf */
    for(i=0;i<nbins*NumLos;i++) {
        double tau = exp(-scale*tau_H1[i])*PBINS;
        if(tau >=PBINS-1)
            flux_pdf[PBINS-1]++;
        else
            flux_pdf[(int)floor(tau)]++;
    }
    /*Normalise*/
    for(i=0;i<PBINS;i++){
        flux_pdf[i]/=(NumLos*nbins);
        flux_pdf[i]*=(PBINS);
    }
    return;
}

/*****************************************************************************/
void help()
{
           fprintf(stderr, "Usage: ./statistic -i input spectra -o output_directory\n"
                           "-r turns off rescaling -s statistic: bit 1 is pdf, 2 is power,\n"
                           "3 is transverse power, 4 is bispectrum\n");
           return;
}


