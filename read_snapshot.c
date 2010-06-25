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
#include "global_vars.h"
#include "parameters.h"

/* this routine loads particle data from Gadget's default
 * binary file format. (A snapshot may be distributed
 * into multiple files. */
int load_snapshot(char *fname, int files,int old, pdata *P)
{
  FILE *fd;
  gadget_header *headers;
  int i,k,n,NumRead,ntot_withmasses=0;
  int NumPart[6]={0,0,0,0,0,0};
  /*Loop over the headers first, to get totals, then loop over data*/
  headers=malloc(files*sizeof(struct io_header_1));
  if(!headers)
  {
    fprintf(stderr, "Error allocating memory for headers.\n");
    exit(1);
  }
  /*First read all the headers, allocate some memory and work out the totals.*/
  for(i=0; i<files; i++)
  {
    char buf[200];
    if(files>1)
      sprintf(buf,"%s.%d",fname,i);
    else
      sprintf(buf,"%s",fname);
    
    if(!(fd=fopen(buf,"r")))
     {
   		fprintf(stderr,"Error opening file %s for reading!\n", buf);
   		exit(1);
     }
     if(!read_gadget_head(headers+i, fd, old))
     {
   		fprintf(stderr,"Error reading file header!\n");
   		exit(1);
     }
     //By now we should have the header data.
     fclose(fd);
  }
  atime= headers[0].time;
  redshift= headers[0].redshift;
  box100 = headers[0].BoxSize;
  omega0 = headers[0].Omega0;
  omegaL = headers[0].OmegaLambda;
  h100 = headers[0].HubbleParam;
  /*Check that the headers all match*/
  for(i=0; i<files; i++)
  {
    /*Load in the number of different particle types
     * NumPart[k] holds the total number of particles of that type*/
    for(k=0; k<6; k++)
        NumPart[k] += headers[i].npart[k]; /* Total particle number */
    if(!NumPart[PARTTYPE])
      {
        fprintf(stderr, "No particles of type  %d! Exiting.\n\n",PARTTYPE);
        exit(1);
      }
    for(k=0; k<6; k++)
        if(headers[i].mass[k]==0 && headers[i].npart[k])
          ntot_withmasses+= headers[i].npart[k];
    if(atime != headers[i].time ||
       redshift!= headers[i].redshift ||
       box100 != headers[i].BoxSize ||
       omega0 != headers[i].Omega0 ||
       omegaL != headers[i].OmegaLambda ||
       h100 != headers[i].HubbleParam  || 
       headers[i].mass[1] !=headers[0].mass[1] ) 
       {
            fprintf(stderr, "Error. File headers are inconsistent.\n");
            exit(4);
       }
  }
  printf("NumPart=[%d,%d,%d,%d,%d,%d], ",NumPart[0],NumPart[1],NumPart[2],NumPart[3],NumPart[4],NumPart[5]);
  printf("Masses=[%g %g %g %g %g %g], ",headers[0].mass[0],headers[0].mass[1],headers[0].mass[2],headers[0].mass[3],headers[0].mass[4],headers[0].mass[5]);
  printf("Particles with mass = %d\n\n",ntot_withmasses);
  printf("Redshift=%g, Ω_M=%g Ω_L=%g\n",redshift,omega0,omegaL);
  printf("Expansion factor = %f\n",atime);
  printf("Hubble = %g Box=%g \n",h100,box100);

  if(!(alloc_parts(P,NumPart[PARTTYPE])))
  {
    fprintf(stderr,"failed to allocate memory.\n\n");
    exit(1);
  }
  for(i=0, NumRead=0; i<files; i++)
  {
    char buf[200];
    int Nstart=0;
    int Ntype=0;
    /*Particles *after* the ones we are interested in, for format 1 files*/
    int Nfinish=0;
    if(files>1)
      sprintf(buf,"%s.%d",fname,i);
    else
      sprintf(buf,"%s",fname);
    
    printf("\nReading file : %s\n", buf);
    if(!(fd=fopen(buf,"r")))
    {
      fprintf(stderr,"can't open file `%s`\n\n",buf);
      exit(2);
    }
    /*Seek past the header*/
    if(old) fseek(fd,2*sizeof(int)+sizeof(struct io_header_1),SEEK_CUR);
    /* Now read the particle data */ 
    Ntype = headers[i].npart[PARTTYPE];
    for(n=0; n<PARTTYPE; n++)
            Nstart+=headers[i].npart[n];
    for(n=PARTTYPE+1; n<6; n++)
            Nfinish+=headers[i].npart[n];
    /*Particle positions */
    printf("Reading position and velocity...\n");
    /*Seek past particles we care not about*/
    if(old) fseek(fd,3*Nstart*sizeof(float),SEEK_CUR);
    read_gadget_float3((*P).Pos+3*NumRead,"POS ",Nstart,Ntype,fd,old);
    if(old) fseek(fd,3*(Nfinish+Nstart)*sizeof(float),SEEK_CUR);
    /* Peculiar velocites */
    read_gadget_float3((*P).Vel+3*NumRead,"VEL ",Nstart,Ntype,fd,old);
    /*Seek past the IDs*/
    if(old) fseek(fd,(2+Ntype+Nfinish+Nstart)*sizeof(int)+(3*Nfinish)*sizeof(float),SEEK_CUR);
    /* Particles masses  */
    printf("Reading mass and temp...\n");

    if(headers[i].mass[PARTTYPE])
        for(n=0;n<Ntype;n++)
          (*P).Mass[NumRead+n] = headers[i].mass[PARTTYPE];
    else
        read_gadget_float((*P).Mass+NumRead,"MASS",Nstart,Ntype,fd,old);
    if(headers[i].mass[0])
          omegab = headers[0].mass[PARTTYPE]/(headers[0].mass[PARTTYPE]+headers[0].mass[1])*omega0;
    else
          omegab = (*P).Mass[0]/((*P).Mass[0]+headers[0].mass[1])*omega0;
    /*Seek past the last masses*/
    if(old) fseek(fd,headers[i].npart[4]*sizeof(float),SEEK_CUR);
    if(PARTTYPE == 0)
      { 
        /*The internal energy of all the Sph particles is read in */
        read_gadget_float((*P).U+NumRead,"U   ",Nstart,Ntype,fd,old);
        /*Seek past RHO*/
        if(old) fseek(fd,2*sizeof(int)+Ntype*sizeof(float),SEEK_CUR);
        /* The free electron fraction */
        if(headers[i].flag_cooling)
          {
            printf("Reading electron fractions...\n");
        #ifndef GADGET3
            read_gadget_float((*P).Ne+NumRead,"NHP ",Nstart,Ntype,fd,old);
        #else
            /* Gadget-III changes the block names*/
            read_gadget_float((*P).Ne+NumRead,"NE  ",Nstart,Ntype,fd,old);
        #endif
        #ifdef HELIUM
            read_gadget_float((*P).NHep+NumRead,"NHEP",Nstart,Ntype,fd,old);
        #endif  
            /*Seek past two NHEP blocks*/
            /* The HI fraction, nHI/nH */
            read_gadget_float((*P).NH0+NumRead,"NH  ",Nstart,Ntype,fd,old);
            /*An NHE block*/
          }

        /* The smoothing length */	  
       read_gadget_float((*P).h+NumRead,"HSML",0,Ntype,fd,old);
       /*SFR: Don't need to read this.*/
       /*if(old) fseek(fd,2*sizeof(int)+Ntype*sizeof(float),SEEK_CUR);*/
      }
    NumRead+=Ntype;
    fclose(fd);
  }
  printf("P[%d].Pos = [%g %g %g]\n", 0, (*P).Pos[0], (*P).Pos[1],(*P).Pos[2]);
  printf("P[%d].Vel = [%g %g %g]\n", 0, (*P).Vel[0], (*P).Vel[1],(*P).Vel[2]);
  printf("P[%d].Mass = %e Ω_B=%g\n\n", NumRead, (*P).Mass[NumRead-1],omegab);
  printf("P[%d].U = %f\n\n", NumRead, (*P).U[NumRead-1]);
  printf("P[%d].Ne = %e\n", NumRead, (*P).Ne[NumRead-1]);
  printf("P[%d].NH0 = %e\n", NumRead, (*P).NH0[NumRead-1]);
  printf("P[%d].h = %f\n",NumRead, (*P).h[NumRead-1]);
  free(headers);
  return NumRead;
}

int alloc_parts(pdata* P, int np)
{
    return ((*P).Vel=malloc(np*3*sizeof(float))) &&
    ((*P).Pos=malloc(np*3*sizeof(float))) &&
    ((*P).Mass=malloc(np*sizeof(float))) &&
    ((*P).U=malloc(np*sizeof(float))) &&
    ((*P).NH0=malloc(np*sizeof(float))) &&
    ((*P).Ne=malloc(np*sizeof(float))) &&
#ifdef HELIUM
    ((*P).NHep=malloc(np*sizeof(float))) &&
#endif
    ((*P).h=malloc(np*sizeof(float)));
}

void free_parts(pdata* P)
{
    free((*P).Vel);
    free((*P).Pos);
    free((*P).Mass);
    free((*P).U);
    free((*P).NH0);
    free((*P).Ne);
    free((*P).h);
#ifdef HELIUM
    free((*P).NHep);
#endif
    return;
}
