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
#include "global_vars.h"
#include "parameters.h"

/*Structure and comparison function for sorted list of los*/
typedef struct _sort_los
{
        int orig_index;
        double priax;
        /*This is xx, unless iaxis=1, in which case it is yy*/
} sort_los;

int compare_xx(const void *a, const void *b)
{
  if(((sort_los *) a)->priax < (((sort_los *) b)->priax))
    return -1;

  if(((sort_los *) a)->priax > (((sort_los *) b)->priax))
    return +1;

  return 0;
}

int get_list_of_near_lines(const double xx,const double yy,const double zz,const double hh, const double boxsize,const los *los_table, const int NumLos,sort_los* sort_los_table,int nxx, int *index_nr_lines, double *dr2_lines);

/*****************************************************************************/
/* This function rescales various things and calculates the absorption*/
#ifndef HELIUM
void Compute_Absorption(double * tau_H1, double *rhoker_H,interp * H1, const double Hz, const double h100, const double box100, const double atime, const double omegab)
#else
void Compute_Absorption(double * tau_H1, double *rhoker_H, interp * H1,double * tau_He2,interp * He2, const double Hz, const double h100, const double box100, const double atime, const double omegab)
#endif
{
  const double H0 = 1.0e5/MPC; /* 100kms^-1Mpc^-1 in SI */ 
    /* Critical matter/energy density at z = 0.0 */
  const double rhoc = 3.0 * (H0*h100)*(H0*h100) / (8.0 * M_PI * GRAVITY); /* kgm^-3 */
  /* Mean hydrogen mass density of the Universe */
  const double critH = (rhoc * omegab * XH) / (atime*atime*atime); /* kgm^-3*/
  /* Conversion factors from internal units */
  const double rscale = (KPC*atime)/h100;   /* convert length to m */
  const double vscale = sqrt(atime);        /* convert velocity to kms^-1 */
  const double mscale = (1.0e10*SOLAR_MASS)/h100; /* convert mass to kg */
  const double escale = 1.0e6;           /* convert energy/unit mass to J kg^-1 */
  const double tscale = ((GAMMA-1.0) * HMASS * PROTONMASS * escale ) / BOLTZMANN; /* convert (with mu) T to K */
  /*const double hscale = rscale * 0.5;*/ /* Note the factor of 0.5 for this kernel definition */
  /*    Calculate the length scales to be used in the box */
  const double dzgrid   = (box100) / (double)NBINS; /* bin size (kpc) */
  const double vmax = box100 * Hz * rscale/ MPC; /* box size (kms^-1) */
  const double vmax2 = vmax/2.0; /* kms^-1 */
  const double dvbin = vmax / (double)NBINS; /* bin size (kms^-1) */

  /* Absorption cross-sections m^2 */
  const double sigma_Lya_H1  = sqrt(3.0*M_PI*SIGMA_T/8.0) * LAMBDA_LYA_H1  * FOSC_LYA;
  /* Prefactor for optical depth  */
  const double A_H1 = rscale*sigma_Lya_H1*C*dzgrid/sqrt(M_PI);  
#ifdef HELIUM
  const double sigma_Lya_He2 = sqrt(3.0*M_PI*SIGMA_T/8.0) * LAMBDA_LYA_HE2 * FOSC_LYA;
  const double A_He2 =  sigma_Lya_He2*C*dzgrid/sqrt(M_PI);
#endif
  int i,j;
    for(i = 0;i<NBINS;i++){
      /* If there are no particles in this bin, rhoker will be zero. 
       * In this case, we set temp and veloc arbitrarily to one, 
       * to avoid nans propagating. Zero rho will imply zero absorption 
       * anyway. */
      if((*H1).rho[i]){       
       (*H1).veloc[i]  = vscale*(*H1).veloc[i]/(*H1).rho[i]; /* HI weighted km s^-1 */ 
       (*H1).temp[i]   = tscale*(*H1).temp[i]/(*H1).rho[i]; /* HI weighted K */
       (*H1).rho[i] *= mscale*pow(rscale,-3); /*Put rhoker in m units*/
      }
      else{
        (*H1).veloc[i]=1;
        (*H1).temp[i]=1;
      }
      rhoker_H[i] *= mscale*pow(rscale,-3);
#ifdef HELIUM
      if((*He2).rho[i]){
       (*He2).veloc[i]  = vscale*(*He2).veloc[i]/(*He2).rho[i]; /* HI weighted km s^-1 */ 
       (*He2).temp[i]   = tscale*(*He2).temp[i]/(*He2).rho[i]; /* HI weighted K */
       (*He2).rho[i] *= mscale*pow(rscale,-3); /*Put rhoker in m units*/
      }
      else{
        (*He2).veloc[i]=1;
        (*He2).temp[i]=1;
      }
#endif
    }
    /* Compute the HI Lya spectra */
    for(i=0;i<NBINS;i++){
        for(j=0;j<NBINS;j++)
          {
            double T0,T1,T2,tau_H1j,aa_H1,u_H1,b_H1,profile_H1;
            double vdiff_H1;
            
            u_H1  = dvbin*j*1.0e3;
        #ifdef PECVEL 
            u_H1 +=(*H1).veloc[j]*1.0e3;
        #endif
            /* Note this is indexed with i, above with j! 
             * This is the difference in velocities between two clouds 
             * on the same sightline*/
            vdiff_H1  = fabs(dvbin*i*1.0e3 - u_H1); /* ms^-1 */
         #ifdef PERIODIC  
      	  if (vdiff_H1 > (vmax2*1.0e3))
      	    vdiff_H1 = (vmax*1.0e3) - vdiff_H1;
         #endif
            b_H1   = sqrt(2.0*BOLTZMANN*(*H1).temp[j]/(HMASS*PROTONMASS));
            T0 = pow(vdiff_H1/b_H1,2);
            T1 = exp(-T0);
            /* Voigt profile: Tepper-Garcia, 2006, MNRAS, 369, 2025 */ 
          #ifdef VOIGT
            aa_H1 = GAMMA_LYA_H1*LAMBDA_LYA_H1/(4.0*M_PI*b_H1);
            T2 = 1.5/T0;	
            if(T0 < 1.0e-6)
              profile_H1  = T1;
            else
              profile_H1  = T1 - aa_H1/sqrt(M_PI)/T0 
                *(T1*T1*(4.0*T0*T0 + 7.0*T0 + 4.0 + T2) - T2 -1.0);
          #else   
            profile_H1 = T1;
          #endif
            tau_H1j  = A_H1  * (*H1).rho[j]  * profile_H1 /(HMASS*PROTONMASS*b_H1) ;
            tau_H1[i]  += tau_H1j;
          }
    }             /* Spectrum convolution */
    /* Compute the HeI Lya spectra: Probably doesn't work now */
#ifdef HELIUM
      for(i=0;i<NBINS;i++)
	{
	  for(j=0;j<NBINS;j++)
	    {
              double T3,T4,T5,tau_He2j,aa_He2,u_He2,b_He2,profile_He2;
              double vdiff_He2;
	      
              /* Note this is indexed with i, above with j! 
               * This is the difference in velocities between two clouds 
               * on the same sightline*/
              u_He2 = dvbin*j*1.0e3;
           #ifdef PECVEL
              u_He2 += (*He2).veloc[j]*1.0e3;
           #endif
              vdiff_He2 = fabs(dvbin*i*1.0e3 - u_He2); /* ms^-1 */
	      
           #ifdef PERIODIC  
		  if (vdiff_He2 > (vmax2*1.0e3))
		     vdiff_He2 = (vmax*1.0e3) - vdiff_He2; 
           #endif
	      
	      b_He2  = sqrt(2.0*BOLTZMANN*(*He2).temp[j]/(HEMASS*PROTONMASS));
	      T3 = (vdiff_He2/b_He2)*(vdiff_He2/b_He2);
	      T4 = exp(-T3);
	      
	      /* Voigt profile: Tepper-Garcia, 2006, MNRAS, 369, 2025 */ 
#ifdef VOIGT
		  aa_He2 = GAMMA_LYA_HE2*LAMBDA_LYA_HE2/(4.0*M_PI*b_He2);
		  T5 = 1.5/T3; 
		   if(T3 < 1.0e-6)
		       profile_He2  = T4;
 		  else
		    profile_He2 = T4 - aa_He2/sqrt(M_PI)/T3 
		    *(T4*T4*(4.0*T3*T3 + 7.0*T3 + 4.0 + T5) - T5 -1.0);
#else   
		  profile_He2 = T4;
#endif
	      tau_He2j = A_He2 * (*He2).rho[j] * profile_He2 /(HMASS*PROTONMASS*b_He2);
	      tau_He2[i] += tau_He2j;
	    }
	}             /* HeI Spectrum convolution */
#endif //HELIUM
      
   for(i = 0;i<NBINS;i++)
   {
      (*H1).rho[i]      /= rhoker_H[i];  /* HI/H */
   #ifdef HELIUM
      (*He2).rho[i]      /= rhoker_H[i];  /* HI/H */
   #endif /*HELIUM*/
      rhoker_H[i]     = log10(mscale*rhoker_H[i]/critH);   /* log H density normalised by mean 
                                                       H density of universe */
   }
  return;
}

/*The size of the thread cache to use below*/
#define CACHESZ 128

/*****************************************************************************/
/*This function does the hard work of looping over all the particles*/
#ifndef HELIUM
void SPH_Interpolation(double * rhoker_H, interp * H1, const int Particles, const int NumLos,const double boxsize, const los *los_table, const pdata *P)
#else
void SPH_Interpolation(double * rhoker_H, interp * H1, interp * He2, const int Particles, const int NumLos,const double boxsize, const los *los_table, const pdata *P)
#endif
{
      /* Loop over particles in LOS and do the SPH interpolation */
      /* This first finds which particles are near this sight line. 
       * Probably a faster way to do that. 
       * Then adds the total density, temp. and velocity for near particles to 
       * the binned totals for that sightline*/
    const double zmingrid = 0.0;
    const double dzgrid   = (boxsize-zmingrid) / (double)NBINS; /* bin size (kpc) */
    const double dzinv    = 1. / dzgrid;
    const double box2     = 0.5 * boxsize;
    int i;
   /*Make a table with a bit more indirection, so we can sort it*/
   sort_los sort_los_table[NumLos];
   /*Need a pointer to the separate structure for los with iaxis=1*/
   sort_los *sort_los_table_xx;
   int nxx=0,nother=0;
   for(i=0;i<NumLos;i++){
       if(los_table[i].axis==1){
             sort_los_table[NumLos-1-nxx].orig_index=i;
             sort_los_table[NumLos-1-nxx].priax=los_table[i].yy;
             nxx++;
       }else{
             sort_los_table[nother].orig_index=i;
             sort_los_table[nother].priax=los_table[i].xx;
             nother++;
       }
   }
   sort_los_table_xx=sort_los_table+NumLos-nxx;
   /*Sort the tables: now the table is sorted we can use bsearch to find the element we are looking for*/
   qsort(sort_los_table,NumLos-nxx,sizeof(sort_los),compare_xx);
   qsort(sort_los_table_xx,nxx,sizeof(sort_los),compare_xx);

  #pragma omp parallel
  {
  /*This is a small thread-local caching mechanism, to avoid
   * massive deadlock around the omp critial section.*/
   int cindex=0;
   double rho_H1[CACHESZ]={0}, temp_H1[CACHESZ]={0},veloc_H1[CACHESZ]={0};
#ifdef HELIUM
   double rho_He2[CACHESZ]={0}, temp_He2[CACHESZ]={0},veloc_He2[CACHESZ]={0};
#endif
   double rho_H[CACHESZ]={0};
   int bins[CACHESZ];
    #pragma omp for
    for(i=0;i<Particles;i++)
    {
      int ind;
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
      int index_nr_lines[NumLos];
      double dr2_lines[NumLos];
      int num_nr_lines=get_list_of_near_lines(xx,yy,zz,hh,boxsize,los_table,NumLos,sort_los_table,nxx,index_nr_lines,dr2_lines);
      for(ind=0;ind<num_nr_lines;ind++)
      {
          int iproc=index_nr_lines[ind];
          double dr2=dr2_lines[ind];
          int iz,ioff,j,iiz;
          /*Load a sightline from the table.*/
          const int iaxis = los_table[iproc].axis;
          const double H1frac = (*P).NH0[i]; /* nHI/nH */
#ifdef HELIUM
          const double He2frac = (*P).NHep[i]; /* nHeII/nH */
#endif
          const double hinv2 = 1. / h2; /* 1/h^2 */
          const double hinv3 = hinv2 / hh; /* 1/h^3 */
          
          const double vr = (*P).Vel[3*i+iaxis-1]; /* peculiar velocity in GII units */
          const double mu = 1.0/(XH*(0.75+(*P).Ne[i]) + 0.25);
          const double temp = (*P).U[i]*mu; /* T in some strange units */
          double dzmax,zgrid;
	     
	     /* Central vertex to contribute to */
	     if (iaxis == 1)
	       iz = (xx - zmingrid) * dzinv +1  ;
	     else if (iaxis == 2) 
	       iz = (yy - zmingrid) * dzinv +1 ;
	     else 
	       iz = (zz - zmingrid) * dzinv +1;
	     
	     dzmax = sqrt(fabs(h4 - dr2));
	     ioff = (int)(dzmax * dzinv) +1;
	     
	     /* Loop over contributing vertices */
	     for(iiz = iz-ioff; iiz < iz+ioff+1 ; iiz++)
	       {
                 double deltaz,dz,q,kernel,velker,temker;
	         j = iiz;
	         j = ((j-1+10*NBINS) % NBINS);
	         
	         zgrid = zmingrid + (double)(j) * dzgrid;
	         
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
	        
	        dr2 = dr2 + (dz*dz);
	        if (dr2 > h4)
	  	   continue;
	        q = sqrt(dr2 * hinv2);
	        if (q <= 1.)
	          kernel = (1.+ (q*q) * (-1.5 + 0.75 * q) )/M_PI;
	        else
	          kernel = 0.25*(2.0-q)*(2.0-q)*(2.0-q)/M_PI;
	        
	        kernel *= hinv3; 

	        kernel *= (*P).Mass[i]; /* kg (kpc)^-3 */
	        velker = vr * kernel; /* kg (kpc)^-3 * km s^-1 */
	        temker = temp * kernel; /* kg (kpc)^-3 * K */
                /*Only one thread can update the global memory at a time.
                 * This adds a small thread-local cache.
                 * Add stuff to the cache*/
                bins[cindex]=iproc*NBINS+j;
                rho_H[cindex]  += kernel * XH;
	        rho_H1[cindex] += kernel * XH * H1frac;
	        veloc_H1[cindex] += velker * XH * H1frac;
	        temp_H1[cindex] += temker * XH * H1frac;
                #ifdef HELIUM
                  rho_He2[cindex] += kernel * XH * He2frac;
                  veloc_He2[cindex] += velker * XH * He2frac;
                  temp_He2[cindex] += temker * XH * He2frac;
                #endif
                cindex++;
                /*Empty the cache when it is full
                 * This is a critical section*/
                if(cindex == CACHESZ){
                  #pragma omp critical
                  {
                          for(cindex=0;cindex<CACHESZ;cindex++){
                             rhoker_H[bins[cindex]]  += rho_H[cindex];
	                     (*H1).rho[bins[cindex]] += rho_H1[cindex];
	                     (*H1).veloc[bins[cindex]] += veloc_H1[cindex];
	                     (*H1).temp[bins[cindex]] +=temp_H1[cindex];
                           #ifdef HELIUM
	                     (*He2).rho[bins[cindex]] += rho_He2[cindex];
	                     (*He2).veloc[bins[cindex]] += veloc_He2[cindex];
	                     (*He2).temp[bins[cindex]] +=temp_He2[cindex];
                           #endif
                          }
                  }/*End critical block*/
                  /* Zero the cache at the end*/
                  for(cindex=0;cindex<CACHESZ;cindex++){
                          rho_H[cindex] = rho_H1[cindex] = veloc_H1[cindex] = temp_H1[cindex] = 0;
                  #ifdef HELIUM
                          rho_He2[cindex] = veloc_He2[cindex] = temp_He2[cindex] = 0;
                  #endif
                  }
                  cindex=0;
                }
	       }        /* loop over contributing vertices */
	}  /*Loop over LOS*/               
    } /* Loop over particles*/
    /*Also empty the cache at the end*/
    #pragma omp critical
    {
            for(i=0;i<cindex;i++){
               rhoker_H[bins[i]]  += rho_H[i];
               (*H1).rho[bins[i]] += rho_H1[i];
               (*H1).veloc[bins[i]] += veloc_H1[i];
               (*H1).temp[bins[i]] +=temp_H1[i];
             #ifdef HELIUM
               (*He2).rho[bins[i]] += rho_He2[i];
               (*He2).veloc[bins[i]] += veloc_He2[i];
               (*He2).temp[bins[i]] +=temp_He2[i];
             #endif
            }
    }/*End critical block*/
  }/*End parallel*/
    return;
}

/*This implements binary search for a given xx.
 * Returns the index of the last element where xx >= priax */
int find_index(double xx, sort_los* sort_los_table, const int NumLos)
{
        int low,high,mid;
        low=0;
        high=NumLos-1;
        while(high - low > 1)
        {
            mid = (high + low) / 2;
            if(xx < sort_los_table[mid].priax)
                high = mid;
            else
                low = mid;
        }
        return low;
}

int get_near_lines_2nd_axis(const double xx,const double yy,const double zz,const double h4, const double boxsize,const sort_los *sort_los_table, const los *los_table, int *index_nr_lines, double *dr2_lines, const int low, const int high)
{
      int ind,num_nr_lines=0;
      for(ind=low;ind<high;ind++)
      {
          double dr,dr2;
          const int iproc=sort_los_table[ind].orig_index;
          /*Load a sightline from the table.*/
          const int iaxis = los_table[iproc].axis;
          double xproj,yproj,zproj;

          xproj = los_table[iproc].xx;
          yproj = los_table[iproc].yy;
          zproj = los_table[iproc].zz;

	  /*    Distance to projection axis */
	  if (iaxis == 1)
	    dr = fabs(yy-yproj);
          else
	    dr = fabs(xx-xproj);

          dr2 = dr*dr;

	  if (iaxis == 3)
            dr = fabs(yy - yproj);
          else
	    dr = fabs(zz - zproj);

	  if (dr > 0.5*boxsize)
	    dr = boxsize - dr; /* between 0 and box/2 */

	  dr2 = dr2 + (dr*dr);

          /*If close in the second coord, save line*/
	  if (dr2 <= h4){
                  index_nr_lines[num_nr_lines]=iproc;
                  dr2_lines[num_nr_lines]=dr2;
                  num_nr_lines++;
          }
      }
      return num_nr_lines;
}

/*This function takes a particle position and returns a list of the indices of lines near it in index_nr_lines
 * Near is defined as: dx^2+dy^2 < 4h^2 */
int get_list_of_near_lines(const double xx,const double yy,const double zz,const double hh, const double boxsize,const los *los_table, const int NumLos,sort_los* sort_los_table,int nxx, int *index_nr_lines, double *dr2_lines)
{
      const double h4 = 4.*hh*hh;           /* 2 smoothing lengths squared */
      int low,high;
      int num_nr_lines=0;
      double ff;
      /*Need a pointer to the separate structure for los with iaxis=1*/
      sort_los *sort_los_table_xx;
      sort_los_table_xx=sort_los_table+NumLos-nxx;
      if(nxx < NumLos){
        /*Now find the elements where dr < 2 hh, wrapping with respect to boxsize*/
        /* First find highest index where xx + 2 hh > priax */
        ff=xx+2*hh;
        if(ff > boxsize)
                ff-=boxsize;
        high=find_index(ff,sort_los_table,NumLos-nxx);
        /* Now find lowest index in what remains where xx - 2 hh < priax */
        ff=xx-2*hh;
        if(ff < 0)
                ff+=boxsize;
        low=find_index(ff,sort_los_table,NumLos-nxx);
        /*This should be the case unless wrapping has occurred*/
        if(low < high)
          num_nr_lines+=get_near_lines_2nd_axis(xx,yy,zz,h4, boxsize,sort_los_table, los_table, index_nr_lines+num_nr_lines, dr2_lines, low, high);
        else{
          num_nr_lines+=get_near_lines_2nd_axis(xx,yy,zz,h4, boxsize,sort_los_table, los_table, index_nr_lines+num_nr_lines, dr2_lines, 0, high);
          num_nr_lines+=get_near_lines_2nd_axis(xx,yy,zz,h4, boxsize,sort_los_table, los_table, index_nr_lines+num_nr_lines, dr2_lines, low, NumLos-nxx);
        }
      }
      if(nxx > 0){
        /*Do the same thing with the table where iaxis=1*/
        /*Now find the elements where dr < 2 hh, wrapping with respect to boxsize*/
        /* First find highest index where xx + 2 hh > priax */
        ff=yy+2*hh;
        if(ff > boxsize)
                ff-=boxsize;
        high=find_index(ff,sort_los_table_xx,nxx);
        /* Now find lowest index in what remains where xx - 2 hh < priax */
        ff=yy-2*hh;
        if(ff < 0)
                ff+=boxsize;
        low=find_index(ff,sort_los_table_xx,nxx);
        /*This should be the case unless wrapping has occurred*/
        if(low < high)
          num_nr_lines+=get_near_lines_2nd_axis(xx,yy,zz,h4, boxsize,sort_los_table_xx, los_table, index_nr_lines+num_nr_lines, dr2_lines, low, high);
        else{
          num_nr_lines+=get_near_lines_2nd_axis(xx,yy,zz,h4, boxsize,sort_los_table_xx, los_table, index_nr_lines+num_nr_lines, dr2_lines, 0, high);
          num_nr_lines+=get_near_lines_2nd_axis(xx,yy,zz,h4, boxsize,sort_los_table_xx, los_table, index_nr_lines+num_nr_lines, dr2_lines, low, nxx);
        }
      }
      return num_nr_lines;
}
