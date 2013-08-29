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
#include "index_table.h"
#include "gadgetreader.hpp"


/* Snapshot information */
#define PARTTYPE 0/* Particle type required */

int load_header(const char *fname,double  *atime, double *redshift, double * Hz, double *box100, double *h100)
{
#ifdef GADGET3
  GadgetReader::GSnap snap(fname);
#else
  const std::string blocks[14]={"HEAD","POS ","VEL ","ID  ","MASS","U   ","RHO ","NHP ","NHEP","NHEQ","NH  ","NHE ","HSML","SFR "};
  std::vector<std::string> BlockNames(blocks,blocks+14);
  GadgetReader::GSnap snap(fname, true,&BlockNames);
#endif
  if(snap.GetNumFiles() <= 0)
          return -1;
  (*atime)= snap.GetHeader().time;
  (*redshift)= snap.GetHeader().redshift;
  (*box100) = snap.GetHeader().BoxSize;
  (*h100) = snap.GetHeader().HubbleParam;
  (*Hz)=100.0*(*h100) * sqrt(1.+snap.GetHeader().Omega0*(1./(*atime)-1.)+snap.GetHeader().OmegaLambda*((pow(*atime,2)) -1.))/(*atime);
  printf("NumPart=[%ld,%ld,%ld,%ld,%ld,%ld), ",snap.GetNpart(0),snap.GetNpart(1),snap.GetNpart(2),snap.GetNpart(3),snap.GetNpart(4),snap.GetNpart(5));
  printf("Masses=[%g %g %g %g %g %g], ",snap.GetHeader().mass[0],snap.GetHeader().mass[1],snap.GetHeader().mass[2],snap.GetHeader().mass[3],snap.GetHeader().mass[4],snap.GetHeader().mass[5]);
  printf("Redshift=%g, Ω_M=%g Ω_L=%g\n",(*redshift),snap.GetHeader().Omega0,snap.GetHeader().OmegaLambda);
  printf("Expansion factor = %f\n",(*atime));
  printf("Hubble = %g Box=%g \n",(*h100),(*box100));
  return 0;


}

/* this routine loads particle data from Gadget's default
 * binary file format. (A snapshot may be distributed
 * into multiple files. */
int64_t load_snapshot(const char *fname,int64_t StartPart,int64_t MaxRead, pdata *P, double *omegab)
{
#ifdef GADGET3
  GadgetReader::GSnap snap(fname);
#else
  const std::string blocks[14]={"HEAD","POS ","VEL ","ID  ","MASS","U   ","RHO ","NHP ","NHEP","NHEQ","NH  ","NHE ","HSML","SFR "};
  std::vector<std::string> BlockNames(blocks,blocks+14);
  GadgetReader::GSnap snap(fname, true,&BlockNames);
#endif
  int64_t NumPart;
  if(MaxRead > 0)
        NumPart = std::min(MaxRead,snap.GetNpart(PARTTYPE)-StartPart);
  else
        NumPart = snap.GetNpart(PARTTYPE)-StartPart;
  if(NumPart ==0)
          return 0;
  if(!(alloc_parts(P,NumPart)))
  {
    fprintf(stderr,"failed to allocate memory.\n\n");
    exit(1);
  }
  printf("Reading from %ld to %ld\n",StartPart,StartPart+NumPart);
  snap.GetBlock("POS ",(*P).Pos,NumPart,StartPart, (1<<N_TYPE)-1-(1<<PARTTYPE));
  snap.GetBlock("VEL ",(*P).Vel,NumPart,StartPart, (1<<N_TYPE)-1-(1<<PARTTYPE));
  /* Particles masses  */
  if(snap.GetHeader().mass[PARTTYPE])
        for(int i=0; i< NumPart;i++)
           (*P).Mass[i] = snap.GetHeader().mass[PARTTYPE];
  else{
        /*Set up types to skip; skip all types not in the mass array*/
        int skip = (1<<N_TYPE)-1-(1<<PARTTYPE);
        for(int i=0; i< N_TYPE; i++)
           if(snap.GetHeader().mass[i])
                   skip-=1<<i;
        snap.GetBlock("MASS",(*P).Mass,NumPart,StartPart, skip);
  }
  for(int i=0; i< NumPart;i++)
  if ((*P).Mass[i] != (*P).Mass[0]){
        fprintf(stderr, "i=%d N = %ld Mass change: %g\n",i,NumPart, (*P).Mass[i]);
        exit(0);
  }
  (*omegab) = (*P).Mass[0]/((*P).Mass[0]+snap.GetHeader().mass[1])*snap.GetHeader().Omega0;
  /*Seek past the last masses*/
  if(PARTTYPE == 0)
    { 
      /*The internal energy of all the Sph particles is read in */
      snap.GetBlock("U   ",(*P).U,NumPart,StartPart,0);
      /* The free electron fraction */
      if(snap.GetHeader().flag_cooling)
        {
          /* Some versions of Gadget have Ne, some have NHP, NHEP and NHEPP,
           * which I map to NHEQ.*/
          /* Use that the universe is neutral, so 
           * NE = NHP + NHEP +2 NHEPP*/
      #ifndef GADGET3
          int k;
          snap.GetBlock("NHP ",(*P).Ne,NumPart,StartPart,0);
          /*Use the space for HSML as temp space*/
          snap.GetBlock("NHEP",(*P).h,NumPart,StartPart,0);
          for(k=0;k<NumPart;k++){
                  (*P).Ne[k]+=(*P).h[k];
          }
          snap.GetBlock("NHEQ",(*P).h,NumPart,StartPart,0);
          for(k=0;k<NumPart;k++){
                  (*P).Ne[k]+=2*(*P).h[k];
          }
      #else
          snap.GetBlock("NE  ",(*P).Ne,NumPart,StartPart,0);
      #endif
      #ifdef HELIUM
          snap.GetBlock("NHE ",(*P).NHep,NumPart,StartPart,0);
      #endif 
          /* The HI fraction, nHI/nH */
          snap.GetBlock("NH  ",(*P).fraction,NumPart,StartPart,0);
          /*An NHE block*/
        }
     /* The smoothing length */
     snap.GetBlock("HSML",(*P).h,NumPart,StartPart,0);
    }

  if(StartPart==0){
        printf("P[%d].Pos = [%g %g %g]\n", 0, (*P).Pos[0], (*P).Pos[1],(*P).Pos[2]);
        printf("P[%d].Vel = [%g %g %g]\n", 0, (*P).Vel[0], (*P).Vel[1],(*P).Vel[2]);
        printf("P[%ld].Mass = %e Ω_B=%g\n\n", NumPart, (*P).Mass[0],(*omegab));
        printf("P[%ld].U = %f\n\n", NumPart, (*P).U[NumPart-1]);
        printf("P[%ld].Ne = %e\n", NumPart, (*P).Ne[NumPart-1]);
        printf("P[%ld].NH0 = %e\n", NumPart, (*P).fraction[NumPart-1]);
        printf("P[%ld].h = %f\n",NumPart, (*P).h[NumPart-1]);
  }
  return NumPart;
}

