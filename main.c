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


/* Here we load a snapshot file. It can be distributed
 * onto several files (for files>1) */
/**********************************************************************/
int main(int argc, char **argv)
{
  int Npart, NumLos=0, files=0;
  FILE *output;
  int rescale=1;
  int old=0;
  los *los_table=NULL;
  char *ext_table=NULL;
  char *outname=NULL;
  char *outdir=NULL;
  char *indir=NULL;
  char c;
  struct particle_data P;
#ifndef RAW_SPECTRA
  double obs_flux,scale=1;
  double flux_power[(NBINS+1)/2];
  int j,jj;
#endif
  /*Make sure stdout is line buffered even when not 
   * printing to a terminal but, eg, perl*/
  setlinebuf(stdout);
  while((c = getopt(argc, argv, "f:o:i:t:n:rh1")) !=-1)
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
        case 't':
           ext_table=optarg;
           break;
        case '1':
           old=1;
           fprintf(stderr, "WARNING: Reading old-style files is not really supported and may not work\n");
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
  los_table=malloc(NumLos*sizeof(los));
  if(!los_table){
          fprintf(stderr, "Error allocating memory for sightline table\n");
          exit(2);
  }
  Npart=load_snapshot(indir, files, old, &P);
  populate_los_table(los_table,NumLos, ext_table, box100); 
  InitLOSMemory(NumLos);
  if(!PARTTYPE)
    SPH_interpolation(NumLos,Npart, los_table, &P);
  /*Free the particle list once we don't need it*/
  free_parts(&P);
  /*If the spectra flag is set, output the raw spectra to a file, 
   * instead of the flux power spectrum directly.*/
#ifdef RAW_SPECTRA
  /*Output to a file*/
  if(!(outname=malloc((strlen(outdir)+29)*sizeof(char))) || !strcpy(outname,outdir) || !(outname=strcat(outname, "_spectra.dat")))
  {
    fprintf(stderr, "Some problem with file output strings\n");
    exit(1);
  }
  if(!(output=fopen(outname,"w")))
  {
          fprintf(stderr, "Error opening %s: %s\n",outname, strerror(errno));
  }
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
#if 0
  fwrite(posaxis,sizeof(double),NBINS,output);          /* pixel positions, comoving kpc/h */
  fwrite(velaxis,sizeof(double),NBINS,output);          /* pixel positions, km s^-1 */
#endif
  fclose(output);
#else /*If outputting flux power directly*/
  /*Calculate mean flux*/
  /*Changing mean flux by a factor of ten changes the P_F by a factor of three*/
  /* Do this even if not rescaling, as we want to know the mean flux*/
  obs_flux= exp(-TAU_EFF);
  scale=mean_flux(tau_H1, NBINS*NumLos,obs_flux,0.001 );
  printf("scale=%g\n",scale);
  /*If no rescale, we output the non-rescaled power spectrum as well*/
  if(!rescale){
     calc_power_spectra(flux_power,tau_H1,1.0,0,NumLos);
     printf("Outputting average flux power spectrum\n");
     outname=malloc((strlen(outdir)+25)*sizeof(char));
     if(!strcpy(outname,outdir) || !(outname=strcat(outname, "_no_rescale_flux_power.txt")))
     {
       fprintf(stderr, "Some problem with the strings\n");
       exit(1);
     }
     output=fopen(outname,"w");
     for(j=0; j<(NBINS+1)/2;j++)
     {
             /*First value is k*/
        fprintf(output, "%g %g\n", j+0.5, flux_power[j]);
     }
     fclose(output);
  }
  calc_power_spectra(flux_power,tau_H1,scale,TAU_EFF,NumLos);
  /*Average the power spectrum*/
  printf("Outputting non-rescaled average flux power spectrum\n");
  if(!(outname=malloc((strlen(outdir)+25)*sizeof(char))) || !strcpy(outname,outdir) || !(outname=strcat(outname, "_flux_power.txt")))
  {
    fprintf(stderr, "Some problem with file output strings\n");
    exit(1);
  }
  if(!(output=fopen(outname,"w")))
  {
          fprintf(stderr, "Error opening %s: %s\n",outname, strerror(errno));
  }
  for(j=0; j<(NBINS+1)/2;j++)
  {
          /*First value is k*/
     fprintf(output, "%g %g\n", j+0.5, flux_power[j]);
  }
  fclose(output);
  #endif /* RAW_SPECTRA*/
  /*Free other memory*/
  free(los_table);
  free(outname);
  FreeLOSMemory();
  return 0;
}
/**********************************************************************/

void help()
{
           fprintf(stderr, "Usage: ./extract -f NUMFILES -n NUMLOS -i filename (ie, without the .0) -o output_file (_flux_power.txt or _spectra.dat will be appended)\n"
                  "-t table_file will read line of sight coordinates from a table.\n"
                  "-r turns off mean flux rescaling\n");
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
                        if(axis > 3 || axis <0 || xx > boxm || xx < 0 || 
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

