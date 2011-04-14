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
#include "gadgetreader.hpp"

/* this routine loads particle data from Gadget's default
 * binary file format. (A snapshot may be distributed
 * into multiple files. */
int load_snapshot(char *fname, pdata *P,
  double  *atime, double *redshift, double * Hz, double *box100, double *h100, double *omegab)
{
  gadget_header header;
  GadgetReader::Gsnap snap(fname);
  int64_t NumPart;
  (*atime)= snap.GetHeader().time;
  (*redshift)= snap.GetHeader().redshift;
  (*box100) = snap.GetHeader().BoxSize;
  (*h100) = snap.GetHeader().HubbleParam;
  (*Hz)=100.0*(*h100) * sqrt(1.+snap.GetHeader().Omega0*(1./(*atime)-1.)+headers[0].OmegaLambda*((pow(*atime,2)) -1.))/(*atime);
  NumPart = snap.GetNpart(PARTTYPE);
  if(NumPart ==0)
          return 0;
  printf("NumPart=[%d,%d,%d,%d,%d,%d), ",snap.GetNpart(0),snap.GetNpart(1),snap.GetNpart(2),snap.GetNpart(3),snap.GetNpart(4),snap.GetNpart(5));
  printf("Masses=[%g %g %g %g %g %g], ",snap.GetHeader().mass[0],snap.GetHeader().mass[1],snap.GetHeader().mass[2],snap.GetHeader().mass[3],snap.GetHeader().mass[4],snap.GetHeader().mass[5]);
  printf("Particles with mass = %d\n\n",ntot_withmasses);
  printf("Redshift=%g, Ω_M=%g Ω_L=%g\n",(*redshift),snap.GetHeader().Omega0,snap.GetHeader().OmegaLambda);
  printf("Expansion factor = %f\n",(*atime));
  printf("Hubble = %g Box=%g \n",(*h100),(*box100));

  if(!(alloc_parts(P,NumPart[PARTTYPE])))
  {
    fprintf(stderr,"failed to allocate memory.\n\n");
    exit(1);
  }
  printf("Reading position and velocity...\n");
  snap.GetBlock("POS ",(*P).Pos,NumPart,0, (1<<N_TYPE)-1-(1<<PARTTYPE));
  snap.GetBlock("VEL ",(*P).Vel,NumPart,0, (1<<N_TYPE)-1-(1<<PARTTYPE));
  /* Particles masses  */
  printf("Reading mass and temp...\n");
  /*DANGER: This assumes that we are reading baryons and that masses are constant*/
  if(snap.GetHeader().mass[PARTTYPE])
      for(n=0;n<Ntype;n++)
        (*P).Mass = snap.GetHeader().mass[PARTTYPE];
  else
        snap.GetBlock("MASS",(*P).Mass,1,0,0);
  (*omegab) = (*P).Mass/((*P).Mass+snap.GetHeader().mass[1])*snap.GetHeader().Omega0;
  /*Seek past the last masses*/
  if(PARTTYPE == 0)
    { 
      /*The internal energy of all the Sph particles is read in */
      snap.GetBlock("U   ",(*P).U,NumPart,0,0);      
      /* The free electron fraction */
      if(headers[i].flag_cooling)
        {
          float *temp;
          int k;
          printf("Reading electron fractions...\n");
          /* Some versions of Gadget have Ne, some have NHP, NHEP and NHEPP,
           * which I map to NHEQ.*/
          /* Use that the universe is neutral, so 
           * NE = NHP + NHEP +2 NHEPP*/
      #ifndef GADGET3
          snap.GetBlock("NHP ",(*P).Ne,NumPart,0,0);      
          /*Use the space for HSML as temp space*/
          snap.GetBlock("NHEP",(*P).h,NumPart,0,0);      
          read_gadget_float(temp,"NHEP",Nstart,Ntype,fd,old);
          for(k=0;k<NumPart;k++){
                  (*P).Ne[k]+=(*P).h[k];
          }
          snap.GetBlock("NHEQ",(*P).h,NumPart,0,0);      
          for(k=0;k<NumPart;k++){
                  (*P).Ne[k]+=2*(*P).h[k];
          }
      #else
          snap.GetBlock("NE  ",(*P).Ne,NumPart,0,0);      
      #endif
      #ifdef HELIUM
          snap.GetBlock("NHE ",(*P).NHep,NumPart,0,0);      
      #endif 
          /* The HI fraction, nHI/nH */
          snap.GetBlock("NH  ",(*P).NH0,NumPart,0,0);      
          /*An NHE block*/
        }
     /* The smoothing length */
     snap.GetBlock("HSML",(*P).h,NumPart,0,0);      
    }
  printf("P[%d].Pos = [%g %g %g]\n", 0, (*P).Pos[0], (*P).Pos[1],(*P).Pos[2]);
  printf("P[%d].Vel = [%g %g %g]\n", 0, (*P).Vel[0], (*P).Vel[1],(*P).Vel[2]);
  printf("P[%d].Mass = %e Ω_B=%g\n\n", NumRead, (*P).Mass[NumRead-1],(*omegab));
  printf("P[%d].U = %f\n\n", NumRead, (*P).U[NumRead-1]);
  printf("P[%d].Ne = %e\n", NumRead, (*P).Ne[NumRead-1]);
  printf("P[%d].NH0 = %e\n", NumRead, (*P).NH0[NumRead-1]);
  printf("P[%d].h = %f\n",NumRead, (*P).h[NumRead-1]);
  free(headers);
#if 0 
    int i;
  /*   Convert to SI units from GADGET-3 units */
  #pragma omp for schedule(static, 128)
  for(i=0;i<Ntype;i++)
  {
      double mu;
      int ic;
      for(ic=0;ic<3;ic++)
	{
	  (*P).Pos[3*i+ic] *= rscale; /* m, physical */
	  (*P).Vel[3*i+ic] *= vscale; /* km s^-1, physical */
	}
      
      (*P).h[i] *= hscale;   /* m, physical */
/*      (*P).Mass[i] = (*P).Mass[i] * mscale; *//* kg */
      /*We leave mass in GADGET units, to prevent a floating overflow
       * when we have poor resolution. (*P).Mass[i] only affects rhoker, 
       * so we simply rescale rhoker later.*/ 

      /* Mean molecular weight */
      mu = 1.0/(XH*(0.75+(*P).Ne[i]) + 0.25);
      (*P).U[i] *= ((GAMMA-1.0) * mu * HMASS * PROTONMASS * escale ) / BOLTZMANN; /* K */
  }
  #pragma omp master
  {
    printf("Converted units.\n");
  }
  #pragma omp barrier
#endif
  return NumRead;
}

int alloc_parts(pdata* P, int np)
{
    return ((*P).Vel=malloc(np*3*sizeof(float))) &&
    ((*P).Pos=malloc(np*3*sizeof(float))) &&
/*     ((*P).Mass=malloc(np*sizeof(float))) && */
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
/*     free((*P).Mass); */
    free((*P).U);
    free((*P).NH0);
    free((*P).Ne);
    free((*P).h);
#ifdef HELIUM
    free((*P).NHep);
#endif
    return;
}
