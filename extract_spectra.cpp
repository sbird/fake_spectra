/* Copyright (c) 2005, J. Bolton
 *      Modified 2009 by Simeon Bird <spb41@cam.ac.uk>
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
#include <string.h>
#include "index_table.h"
#include "parameters.h"

/* Function to rescale the units of the density temperature and velocity skewers*/
void Rescale_Units(double * rho, double * veloc, double * temp, const int nbins, const double h100, const double atime)
{
  /* Conversion factors from internal units */
  const double rscale = (KPC*atime)/h100;   /* convert length to m */
  const double vscale = sqrt(atime);        /* convert velocity to kms^-1 */
  const double mscale = (1.0e10*SOLAR_MASS)/h100; /* convert mass to kg */
  /*const double hscale = rscale * 0.5;*/ /* Note the factor of 0.5 for this kernel definition */
  /*    Calculate the length scales to be used in the box */

  for(int i = 0;i<nbins;i++){
    /* If there are no particles in this bin, rhoker will be zero.
     * In this case, we set temp and veloc arbitrarily to one,
     * to avoid nans propagating. Zero rho will imply zero absorption
     * anyway. */
    if(rho[i]){
     veloc[i]  *= (vscale/rho[i]); /* HI weighted km s^-1 */
     temp[i]   /= rho[i]; /* HI weighted K */
     rho[i] *= mscale*pow(rscale,-3); /*Put rhoker in m units*/
    }
    else{
      veloc[i]=1;
      temp[i]=1;
    }
  }
  return;
}

/*Convert densities for this species into fraction of total hydrogen density.
 * Also rescale rhoker_H */
void Convert_Density(double * rhoker_H, double * rho, const double h100, const double atime, const double omegab)
{
  /* Conversion factors from internal units */
  const double rscale = (KPC*atime)/h100;   /* convert length to m */
  const double mscale = (1.0e10*SOLAR_MASS)/h100; /* convert mass to kg */
  /*Cosmological factors*/
  const double H0 = 1.0e5/MPC; /* 100kms^-1Mpc^-1 in SI */ 
  /* Critical matter/energy density at z = 0.0 */
  const double rhoc = 3.0 * (H0*h100)*(H0*h100) / (8.0 * M_PI * GRAVITY); /* kgm^-3 */
  /* Mean hydrogen mass density of the Universe */
  const double critH = (rhoc * omegab * XH) / (atime*atime*atime); /* kgm^-3*/
  for(int i = 0;i<NBINS;i++)
  {
     rhoker_H[i] *= mscale*pow(rscale,-3);
     rho[i]      /= rhoker_H[i];  /* HI/H */
     rhoker_H[i]     = log10(rhoker_H[i]/critH);   /* log H density normalised by mean
                                                      H density of universe */
  }
  return;
}

/*The size of the thread cache to use below*/
#define CACHESZ 128

/*****************************************************************************/
/*This function does the hard work of looping over all the particles.
 * Can handle an arbitrary number of species, in the array 'species'. 
 * Number of species is given in nspecies. 
 * The fractional abundance of the species, Z = n / n(H) should stored be in P.fraction
 * Output is stored in interp *species, and is rho, temp and vel, in internal gadget units. 
 */
void SPH_Interpolation(double * rhoker_H, interp * species, const int nspecies, const int nbins, const int Particles, const int NumLos,const double boxsize, const los *los_table, IndexTable& sort_los_table,const pdata *P)
{
      /* Loop over particles in LOS and do the SPH interpolation */
      /* This first finds which sightlines are near the particle using the sorted los table 
       * Then adds the total density, temp. and velocity for near particles to 
       * the binned totals for that sightline. Is O(N_part)*O(log n_los) */
    const double dzgrid   = boxsize / (double)nbins; /* bin size (kpc) */
    const double dzinv    = 1. / dzgrid;
    const double box2     = 0.5 * boxsize;
    //Do the T conversion here for convenience - internal energy is messy
    const double escale = 1.0e6;           // convert energy/unit mass to J kg^-1
    /* convert U (J/kg) to T (K) : U = N k T / (γ - 1)
     * T = U (γ-1) μ m_P / k_B
     * where k_B is the Boltzmann constant
     * γ is 5/3, the perfect gas constant
     * m_P is the proton mass
     * μ is 1 / (mean no. molecules per unit atomic weight) calculated in loop.
     */
    const double tscale = ((GAMMA-1.0) * PROTONMASS * escale ) / BOLTZMANN;

  #pragma omp parallel
  {
  /*This is a small thread-local caching mechanism, to avoid
   * massive deadlock around the omp critial section.*/
   int cindex=0;
   double rho[CACHESZ*nspecies], temp[CACHESZ*nspecies], veloc[CACHESZ*nspecies];
   double rho_H[CACHESZ]={0};
   int bins[CACHESZ];
   memset(rho, 0, CACHESZ*nspecies*sizeof(double));
   memset(temp, 0, CACHESZ*nspecies*sizeof(double));
   memset(veloc, 0, CACHESZ*nspecies*sizeof(double));
    #pragma omp for
    for(int i=0;i<Particles;i++)
    {
      /*     Positions (kpc) */
      const double xx = (*P).Pos[3*i+0];
      const double yy = (*P).Pos[3*i+1];
      const double zz = (*P).Pos[3*i+2];
          
      /* Resolution length (kpc) */
      const double hh = (*P).h[i]*0.5; /*Factor of two in this kernel definition*/
      const double h2 = hh*hh;
      const double h4 = 4.*h2;           /* 2 smoothing lengths squared */
/*       if((Particles <20) ||  ((i % (Particles/20)) ==0)) */
/*              printf("Interpolating particle %d.\n",i); */
      std::map<int, double> nearby=sort_los_table.get_near_lines((*P).Pos+3*i,hh);
      for(std::map<int, double>::iterator it = nearby.begin(); it != nearby.end(); ++it)
      {
          int iproc=it->first;
          double dr2=it->second;
          int iz,ioff;
          /*Load a sightline from the table.*/
          const int iaxis = los_table[iproc].axis;
          const double hinv2 = 1. / h2; /* 1/h^2 */
          const double hinv3 = hinv2 / hh; /* 1/h^3 */
          
          const double vr = (*P).Vel[3*i+iaxis-1]; /* peculiar velocity in GII units */
          /*Mean molecular weight:
           * \mu = 1 / molecules per unit atomic weight
           *     = 1 / (X + Y /4 + E)
           *     where E = Ne * X, and Y = (1-X).
           *     Can neglect metals as they are heavy.
           *     Leading contribution is from electrons, which is already included
           *     [+ Z / (12->16)] from metal species
           *     [+ Z/16*4 ] for OIV from electrons. */
          const double mu = 1.0/(XH*(0.75+(*P).Ne[i]) + 0.25);
          const double p_temp = (*P).U[i]*mu * tscale; /* T in K */
          double dzmax,zgrid;
	     
	     /* Central vertex to contribute to */
	     if (iaxis == 1)
	       iz = xx * dzinv;
	     else if (iaxis == 2) 
	       iz = yy * dzinv;
	     else 
	       iz = zz * dzinv;
	     
	     dzmax = sqrt(fabs(h4 - dr2));
	     ioff = (int)(dzmax * dzinv);
	     
	     /* Loop over contributing vertices */
	     for(int iiz = iz-ioff; iiz < iz+ioff+1; iiz++)
	     {
             double deltaz,dz,dist2,velker,temker, kernel;
	         int j = ((iiz+10*nbins) % nbins);
	         
	         zgrid = (double)(j) * dzgrid;
	         
	        if (iaxis == 1)
                  deltaz = zgrid - xx;
	        else if (iaxis == 2)
                  deltaz = zgrid - yy;
	        else
                  deltaz = zgrid - zz;
               
	        if (deltaz > box2)
                deltaz = deltaz - boxsize;
	        if (deltaz < -box2)
                deltaz = deltaz + boxsize;
	        
	        dz = fabs(deltaz);
	        if(dz > box2)
                dz = boxsize - dz;
	        
	        dist2 = dr2 + (dz*dz);
	        if (dist2 > h4)
	  	        continue;

#ifdef SPH_KERNEL
	        const double q = sqrt(dist2 * hinv2);
	        if (q <= 1.)
                    kernel = (1.+ (q*q) * (-1.5 + 0.75 * q) )/M_PI;
	        else
	            kernel = 0.25*(2.0-q)*(2.0-q)*(2.0-q)/M_PI;
	        kernel *= hinv3; 
#else
            //Top-hat kernel for Arepo
            kernel = hinv3 / M_PI /4. *3.;
#endif
	        kernel *= (*P).Mass[i]; /* kg (kpc)^-3 */
	        velker = vr * kernel; /* kg (kpc)^-3 * km s^-1 */
	        temker = p_temp * kernel; /* kg (kpc)^-3 * K */
            /*Only one thread can update the global memory at a time.
             * This adds a small thread-local cache.
             * Add stuff to the cache*/
            bins[cindex]=iproc*nbins+j;
            if(rhoker_H)
                rho_H[cindex]  = kernel;

            for(int l=0;l<nspecies;l++){
                /*GFM_Metals is the total mass in a metal species per unit gas mass.
                 * So use it directly.*/
                rho[cindex*nspecies+l] = kernel * (*P).fraction[nspecies*i+l];
                veloc[cindex*nspecies +l] = velker * (*P).fraction[nspecies*i+l];
                temp[cindex*nspecies+l] = temker * (*P).fraction[nspecies*i+l];
            }
            cindex++;
            /*Empty the cache when it is full
             * This is a critical section*/
            if(cindex == CACHESZ){
              #pragma omp critical
              {
		if(rhoker_H)
                    for(cindex=0;cindex<CACHESZ;cindex++){
                        rhoker_H[bins[cindex]]  += rho_H[cindex];
                    }
                for(cindex=0;cindex<CACHESZ;cindex++){
                    for(int l = 0; l < nspecies; l++){
                        (*species).rho[bins[cindex]*nspecies+l] += rho[cindex*nspecies+l];
                        (*species).veloc[bins[cindex]*nspecies+l] += veloc[cindex*nspecies+l];
                        (*species).temp[bins[cindex]*nspecies+l] +=temp[cindex*nspecies+l];
                        /* Zero the cache at the end*/
                        rho_H[cindex] = rho[cindex*nspecies+l] = veloc[cindex*nspecies+l] = temp[cindex*nspecies+l] = 0;
                    }
                }
              }/*End critical block*/
              cindex=0;
            }
	     }        /* loop over contributing vertices */
	  }  /*Loop over LOS*/               
    } /* Loop over particles*/
    /*Also empty the cache at the end*/
    #pragma omp critical
    {
        for(int i=0;i<cindex;i++){
            for(int l = 0; l < nspecies; l++){
                (*species).rho[bins[i]*nspecies+l] += rho[i*nspecies+l];
                (*species).veloc[bins[i]*nspecies+l] += veloc[i*nspecies+l];
                (*species).temp[bins[i]*nspecies+l] +=temp[i*nspecies+l];
            }
        }
	    if(rhoker_H)
            for(int i=0;i<cindex;i++)
                rhoker_H[bins[i]]  += rho_H[i];
    }/*End critical block*/
  }/*End parallel*/
    return;
}
