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
int load_snapshot(char *fname, int files, int old)
{
  FILE *fd;
  gadget_header *headers;
  float *temp;
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

  if(!(P=malloc(NumPart[PARTTYPE]*sizeof(struct particle_data))))
  {
    fprintf(stderr,"failed to allocate memory.\n\n");
    exit(1);
  }
  for(i=0, NumRead=0; i<files; i++)
  {
    char buf[200];
    int Nstart=0;
    int Ntype=0;
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
    
    /* Now read the particle data */ 
    Ntype = headers[i].npart[PARTTYPE];
    for(n=0; n<PARTTYPE; n++)
            Nstart+=headers[i].npart[n];
    if(!(temp=malloc(3*(floor(Ntype/3)+2)*sizeof(float))))
    {
                    fprintf(stderr, "Failed to allocate temp memory\n");
                    exit(2);
    }
    /*Particle positions */
    printf("Reading position and velocity...\n");
    /*This is a cheat so that peak memory usage is a bit lower*/
    for(k=0; k<3; k++)
    {
        int chunk=floor(Ntype/3);
        int chunke=chunk;
        if(k==2)
           chunke=Ntype-2*chunk;
        read_gadget_float3(temp,"POS ",Nstart+k*chunk,chunke,fd,old);
            for(n=0;n<chunke;n++)
            {
               P[NumRead+n+k*chunk].Pos[0]=temp[3*n];	
               P[NumRead+n+k*chunk].Pos[1]=temp[3*n+1];	
               P[NumRead+n+k*chunk].Pos[2]=temp[3*n+2];	
            }
    }
    
    /* Peculiar velocites */
    for(k=0; k<3; k++)
    {
        int chunk=floor(Ntype/3);
        int chunke=chunk;
        if(k==2)
           chunke=Ntype-2*chunk;
        read_gadget_float3(temp,"VEL ",Nstart+k*chunk,chunke,fd,old);
            for(n=0;n<chunke;n++)
            {
               P[NumRead+n+k*chunk].Vel[0]=temp[3*n];	
               P[NumRead+n+k*chunk].Vel[1]=temp[3*n+1];	
               P[NumRead+n+k*chunk].Vel[2]=temp[3*n+2];	
            }
    }
    free(temp);
    /*Done with triplets, read single valued things*/
    if(!(temp=malloc(Ntype*sizeof(float))))
    {
                    fprintf(stderr, "Failed to allocate temp memory\n");
                    exit(2);
    }
    
    /* Particles masses  */
    printf("Reading mass and temp...\n");
    if(headers[i].mass[PARTTYPE])
        for(n=0;n<Ntype;n++)
          P[NumRead+n].Mass = headers[i].mass[PARTTYPE];
    else
    {
        read_gadget_float(temp,"MASS",Nstart,Ntype,fd, old);
        for(n=0;n<Ntype;n++)
          P[NumRead+n].Mass= temp[n];
    }
    if(headers[i].mass[0])
          omegab = headers[0].mass[PARTTYPE]/(headers[0].mass[PARTTYPE]+headers[0].mass[1])*omega0;
    else
          omegab = P[0].Mass/(P[0].Mass+headers[0].mass[1])*omega0;
    if(PARTTYPE == 0)
      { 
        /*The internal energy of all the Sph particles is read in */
        read_gadget_float(temp,"U   ",Nstart,Ntype,fd, old);
        for(n=0; n<Ntype;n++)
            P[NumRead+n].U=temp[n];
        /* The free electron fraction */
        if(headers[i].flag_cooling)
          {
            printf("Reading electron fractions...\n");
        #ifndef GADGET3
            read_gadget_float(temp,"NHP ",Nstart,Ntype,fd,old);
        #else
            /* Gadget-III changes the block names*/
            read_gadget_float(temp,"NE  ",Nstart,Ntype,fd,old);
        #endif
            for(n=0; n<Ntype;n++)
               P[NumRead+n].Ne=temp[n];

            /* The HI fraction, nHI/nH */
            read_gadget_float(temp,"NH  ",Nstart,Ntype,fd, old);
            for(n=0; n<Ntype;n++)
              P[NumRead+n].NH0=temp[n];
        #ifdef HELIUM
            read_gadget_float(temp,"NHEP",Nstart,Ntype,fd, old);
            for(n=0; n<Ntype;n++)
               P[NumRead+n].NHep=temp[n];
        #endif  
          }

        /* The smoothing length */	  
       read_gadget_float(temp,"HSML",Nstart,Ntype,fd, old);
       for(n=0; n<Ntype;n++)
          P[NumRead+n].h=temp[n];
        
      }
    NumRead+=Ntype;
    fclose(fd);
    free(temp);
  }
  printf("P[%d].Pos = [%g %g %g]\n", 0, P[0].Pos[0], P[0].Pos[1],P[0].Pos[2]);
  printf("P[%d].Vel = [%g %g %g]\n", 0, P[0].Vel[0], P[0].Vel[1],P[0].Vel[2]);
  printf("P[%d].Mass = %e Ω_B=%g\n\n", NumRead, P[1].Mass,omegab);
  printf("P[%d].U = %f\n\n", NumRead, P[1].U);
  printf("P[%d].Ne = %e\n", NumRead, P[1].Ne);
  printf("P[%d].NH0 = %e\n", NumRead, P[1].NH0);
  printf("P[%d].h = %f\n",NumRead, P[1].h);
  free(headers);
  return NumRead;
}
