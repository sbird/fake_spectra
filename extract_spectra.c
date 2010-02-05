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

#include "headers.h"
#include "global_vars.h"
#include "parameters.h"
#define THREAD_ALLOC 10

double *Delta,*posaxis,*velaxis;
double *n_H1,*veloc_H1,*temp_H1,*tau_H1;
float *flux_power;
#ifdef HELIUM
double *n_He2,*veloc_He2,*temp_He2,*tau_He2;
#endif

void InitLOSMemory(int NumLos);
void FreeLOSMemory(void);

/*****************************************************************************/
void SPH_interpolation(int NumLos, int Ntype)
{
  const double Hz=100.0*h100 * sqrt(1.+omega0*(1./atime-1.)+omegaL*((atime*atime) -1.))/atime;
  const double H0 = 1.0e5/MPC; /* 100kms^-1Mpc^-1 in SI */ 
    /* Critical matter/energy density at z = 0.0 */
  const double rhoc = 3.0 * (H0*h100)*(H0*h100) / (8.0 * M_PI * GRAVITY); /* kgm^-3 */
  /* Mean hydrogen mass density of the Universe */
  const double critH = (rhoc * OMEGAB * XH) / (atime*atime*atime); /* kgm^-3*/
  /* Conversion factors from internal units */
  const double rscale = (KPC*atime)/h100;   /* convert length to m */
  const double vscale = sqrt(atime);        /* convert velocity to kms^-1 */
  const double mscale = (1.0e10*SOLAR_MASS)/h100; /* convert mass to kg */
  const double escale = 1.0e6;           /* convert energy/unit mass to J kg^-1 */
  const double hscale = rscale * 0.5; /* Note the factor of 0.5 for this kernel definition */
  /*    Calculate the length scales to be used in the box */
  const double zmingrid = 0.0;
  const double zmaxgrid = box100*rscale;  /* box sizes in physical m */
  const double dzgrid   = (zmaxgrid-zmingrid) / (double)NBINS; /* bin size (physical m) */
  const double dzinv    = 1. / dzgrid;
  const double boxsize  = zmaxgrid;   
  const double box2     = 0.5 * boxsize;
  const double dzbin = box100/ (double)NBINS; /* bin size (comoving kpc/h) */
  const double vmax = box100 * Hz * rscale/ MPC; /* box size (kms^-1) */
  const double vmax2 = vmax/2.0; /* kms^-1 */
  const double dvbin = vmax / (double)NBINS; /* bin size (kms^-1) */
  const double ztime = 1.0/atime - 1.0;

  /* Absorption cross-sections m^2 */
  const double sigma_Lya_H1  = sqrt(3.0*M_PI*SIGMA_T/8.0) * LAMBDA_LYA_H1  * FOSC_LYA;
  /* Prefactor for optical depth  */
  const double A_H1 = sigma_Lya_H1*C*dzgrid/sqrt(M_PI);  
#ifdef HELIUM
  const double sigma_Lya_He2 = sqrt(3.0*M_PI*SIGMA_T/8.0) * LAMBDA_LYA_HE2 * FOSC_LYA;
  const double A_He2 =  sigma_Lya_He2*C*dzgrid/sqrt(M_PI);
#endif
  int iproc;
  FILE *output;

  InitLOSMemory(NumLos);
  
  srand48(241008); /* random seed generator */
  /*   Initialise distance coordinate for iaxis */
  posaxis[0]=0.0;
  velaxis[0]=0.0;
  
  for(iproc=0;iproc<NBINS-1;iproc++)
    {
      posaxis[iproc+1] = posaxis[iproc] + dzbin; /* comoving kpc/h */
      velaxis[iproc+1] = velaxis[iproc] + dvbin; /* physical km s^-1 */
    }
  
#pragma omp parallel
  { 
    int i;
  /*   Convert to SI units from GADGET-3 units */
  #pragma omp for schedule(static, 128)
  for(i=0;i<Ntype;i++)
  {
      double mu;
      int ic;
      for(ic=0;ic<3;ic++)
	{
	  P[i].Pos[ic] *= rscale; /* m, physical */
	  P[i].Vel[ic] *= vscale; /* km s^-1, physical */
	}
      
      P[i].h *= hscale;   /* m, physical */
      P[i].Mass = P[i].Mass * mscale;   /* kg */

      /* Mean molecular weight */
      mu = 1.0/(XH*(0.75+P[i].Ne) + 0.25);
      P[i].U *= ((GAMMA-1.0) * mu * HMASS * PROTONMASS * escale ) / BOLTZMANN; /* K */
  }
  #pragma omp master
  {
    printf("Converting units...done\n");
  }
  #pragma omp barrier
  
  /*    Generate random coordinates for a point in the box */
  #pragma omp for schedule(static, THREAD_ALLOC)
  for(iproc=0;iproc<NumLos;iproc++)
    { 
      double xproj,yproj,zproj;
      int iaxis,iz,ioff,j,iiz,ii;
      double rhoker_H[NBINS],rhoker_H1[NBINS];
      double velker_H1[NBINS],temker_H1[NBINS];
      double temp_H1_local[NBINS],veloc_H1_local[NBINS], tau_H1_local[NBINS];
      float flux_H1_local[NBINS];
      float flux_power_local[(NBINS+1)/2];
#ifdef HELIUM
      double rhoker_He2[NBINS];
      double velker_He2[NBINS],temker_He2[NBINS];
      double temp_He2_local[NBINS],veloc_He2_local[NBINS], tau_He2_local[NBINS];
#endif
      for(i=0; i<NBINS; i++)
      {
         rhoker_H[i]=0;
         rhoker_H1[i]=0;
         velker_H1[i]=0;
         temker_H1[i]=0;
         temp_H1_local[i]=0;
         veloc_H1_local[i]=0;
         tau_H1_local[i]=0;
      }
      /*Pick a random sightline*/
      do	
      	iaxis = (int)(drand48()*4);
      while (iaxis == 0 || iaxis==4); 
      
      xproj = drand48()*box100*rscale;
      yproj = drand48()*box100*rscale;
      zproj = drand48()*box100*rscale;
     if((iproc %10) ==0)
      printf("Interpolating line of sight %d...done\n",iproc);
      
      /* Loop over particles in LOS and do the SPH interpolation */
      /* This first finds which particles are near this sight line. 
       * Probably a faster way to do that. 
       * Then adds the total density, temp. and velocity for near particles to 
       * the binned totals for that sightline*/
      for(i=0;i<Ntype;i++)
	{
	  double xx,yy,zz,hh,h2,h4,dr,dr2;
          double hinv2,hinv3,vr,Temperature,dzmax,H1frac,zgrid;
        #ifdef HELIUM
          double He2frac;
        #endif
	  /*     Positions (physical m) */
	  xx = P[i].Pos[0];
	  yy = P[i].Pos[1];
	  zz = P[i].Pos[2];
              
	  /* Resolution length (physical m) */
	  hh = P[i].h; 
	  h2 = hh*hh; 
	  h4 = 4.*h2;           /* 2 smoothing lengths squared */
	  
	  /*    Distance to projection axis */	  
	  if (iaxis == 1) 
	    dr = fabs(yy-yproj);
          else if (iaxis ==  2)
	    dr = fabs(xx-xproj);
          else 
            dr = fabs(xx-xproj);
	  
	  if (dr > box2) 
	    dr = boxsize - dr; /* Keep dr between 0 and box/2 */
	  
	  if (dr <= 2.*hh) /* dr less than 2 smoothing lengths */
	    {
	      dr2 = dr*dr;
	      
	      if (iaxis == 1)
		dr = fabs(zz - zproj);
              else if (iaxis == 2)
		dr = fabs(zz - zproj);
              else if(iaxis == 3)
		dr = fabs(yy - yproj);
	      
	      if (dr > box2)  
		dr = boxsize - dr; /* between 0 and box/2 */
              
	      dr2 = dr2 + (dr*dr);
	      
	      if (dr2 <= h4)
		{
	           hinv2 = 1. / h2; /* 1/h^2 */
		   hinv3 = hinv2 / hh; /* 1/h^3 */
		   
		   vr = P[i].Vel[iaxis-1]; /* peculiar velocity in km s^-1 */
		   Temperature = P[i].U; /* T in Kelvin */
		   H1frac = P[i].NH0; /* nHI/nH */ 
                #ifdef HELIUM
                   He2frac = P[i].NHep; /* nHeII/nH */
                #endif
		   
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
                       double deltaz,dz,dist2,q,kernel,velker,temker;
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
		      
		      dist2 = dr2 + (dz*dz);		 

		      if (dist2 <= h4)
			{
			  q = sqrt(dist2 * hinv2);
			  if (q <= 1.)
			    kernel = (1.+ (q*q) * (-1.5 + 0.75 * q) )/M_PI;
			  else
			    kernel = 0.25*(2.0-q)*(2.0-q)*(2.0-q)/M_PI;
			  
			  kernel *= hinv3;  

			  kernel *= P[i].Mass; /* kg m^-3 */
			  velker = vr * kernel; /* kg m^-3 * km s^-1 */
			  temker = Temperature * kernel; /* kg m^-3 * K */
			  
			  rhoker_H[j]  += kernel * XH;		 
			  rhoker_H1[j] += kernel * XH * H1frac;
			  velker_H1[j] += velker * XH * H1frac;
			  temker_H1[j] += temker * XH * H1frac;
                        #ifdef HELIUM
                          rhoker_He2[j] += kernel * XH * He2frac;
                          velker_He2[j] += velker * XH * He2frac;
                          temker_He2[j] += temker * XH * He2frac;
                        #endif

			}      /* dist2 < 4h^2 */
		     }        /* loop over contributing vertices */
		 }           /* dx^2+dy^2 < 4h^2 */
	    }               /* dx < 2h */
	}                  /* Loop over particles in LOS */
      
      for(i = 0;i<NBINS;i++)
	{
	  veloc_H1_local[i]  = velker_H1[i]/rhoker_H1[i]; /* HI weighted km s^-1 */ 
	  temp_H1_local[i]   = temker_H1[i]/rhoker_H1[i]; /* HI weighted K */
      	}
      
      /* Compute the HI Lya spectra */
      for(i=0;i<NBINS;i++)
	{
	  for(j=0;j<NBINS;j++)
	    {
              double T0,T1,T2,tau_H1j,aa_H1,u_H1,b_H1,profile_H1;
              double vdiff_H1;
	      
              u_H1  = velaxis[j]*1.0e3;
          #ifdef PECVEL 
              u_H1 +=veloc_H1_local[j]*1.0e3;
          #endif
              /* Note this is indexed with i, above with j! 
               * This is the difference in velocities between two clouds 
               * on the same sightline*/
	      vdiff_H1  = fabs((velaxis[i]*1.0e3) - u_H1); /* ms^-1 */
           #ifdef PERIODIC  
		  if (vdiff_H1 > (vmax2*1.0e3))
		    vdiff_H1 = (vmax*1.0e3) - vdiff_H1;
           #endif
	      b_H1   = sqrt(2.0*BOLTZMANN*temp_H1_local[j]/(HMASS*PROTONMASS));
	      T0 = (vdiff_H1/b_H1)*(vdiff_H1/b_H1);
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
	      tau_H1j  = A_H1  * rhoker_H1[j]  * profile_H1  /(HMASS*PROTONMASS*b_H1);
	      tau_H1_local[i]  += tau_H1j;
	    }
	}             /* Spectrum convolution */
      /* Compute the HeI Lya spectra: Probably doesn't work now */
#ifdef HELIUM
      for(i = 0;i<NBINS;i++)
	{
	  veloc_He2_local[i]  = velker_He2[i]/rhoker_He2[i]; /* HI weighted km s^-1 */ 
	  temp_He2_local[i]   = temker_He2[i]/rhoker_He2[i]; /* HI weighted K */
      	}
      for(i=0;i<NBINS;i++)
	{
	  for(j=0;j<NBINS;j++)
	    {
              double T3,T4,T5,tau_He2j,aa_He2,u_He2,b_He2,profile_He2;
              double vdiff_He2;
	      
              /* Note this is indexed with i, above with j! 
               * This is the difference in velocities between two clouds 
               * on the same sightline*/
              u_He2 = velaxis[j]*1.0e3;
           #ifdef PECVEL
              u_He2 += veloc_He2_local[j]*1.0e3;
           #endif
              vdiff_He2 = fabs((velaxis[i]*1.0e3) - u_He2); /* ms^-1 */
	      
           #ifdef PERIODIC  
		  if (vdiff_He2 > (vmax2*1.0e3))
		     vdiff_He2 = (vmax*1.0e3) - vdiff_He2; 
           #endif
	      
	      b_He2  = sqrt(2.0*BOLTZMANN*temp_He2_local[j]/(HEMASS*PROTONMASS));
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
	      tau_He2j = A_He2 * rhoker_He2[j] * profile_He2 /(HMASS*PROTONMASS*b_He2);
	      tau_He2_local[i] += tau_He2j;
	    }
	}             /* HeI Spectrum convolution */
#endif //HELIUM
      
      for(i = 0;i<NBINS;i++)
	{
	  ii = i + (NBINS*iproc);
	  
	  Delta[ii]     = log10(rhoker_H[i]/critH);   /* log H density normalised by mean 
                                                          H density of universe */
	  n_H1[ii]      = rhoker_H1[i]/rhoker_H[i];  /* HI/H */
          veloc_H1[ii]    = veloc_H1_local[i]; /* HI weighted km s^-1 */ 
	  temp_H1[ii]   = temp_H1_local[i]; /* HI weighted K */
          tau_H1[ii]    = tau_H1_local[i];
        #ifdef HELIUM
	  n_He2[ii]      = rhoker_He2[i]/rhoker_H[i];  /* HI/H */
          veloc_He2[ii]    = veloc_He2_local[i]; /* HI weighted km s^-1 */ 
	  temp_He2[ii]   = temp_He2_local[i]; /* HI weighted K */
          tau_He2[ii]    = tau_He2_local[i];
        #endif
      	}
    }                /* Loop over numlos random LOS */
  }/*End of parallel block*/
  /*Calculate mean flux*/
  double obs_flux= 0.0023*pow(1.0+ztime,3.65);
  double scale=mean_flux(tau_H1, NBINS*NumLos,obs_flux,0.01 );
  int i;
  for(i=0; i<NBINS*NumLos; i++)
  {
    tau_H1[i]*=scale;
  }
#if 0 
#pragma omp_parallel
  {
  #pragma omp for schedule(static, THREAD_ALLOC)
  for(iproc=0;iproc<NumLos;iproc++)
  {
      /* Calculate flux and flux power spectrum */
      for(i=0; i<NBINS; i++)
      {
         flux_H1_local[i]=exp(-tau_H1_local[i]);
      }
      powerspectrum(NBINS, flux_H1_local, flux_power_local);
      /*All non-thread-local memory writing should happen here*/
      /*Write powerspectrum*/
      for(i=0; i<(NBINS+1)/2;i++)
      {
          ii=i+(NBINS+1)/2*iproc;
          flux_power[ii]=flux_power_local[i];
      }
  }/*End loop*/
  }/*End parallel block*/
#endif
  output = fopen("spectra1024.dat","wb");
  fwrite(&ztime,sizeof(double),1,output);
  fwrite(Delta,sizeof(double),NBINS*NumLos,output);     /* gas overdensity */
  fwrite(n_H1,sizeof(double),NBINS*NumLos,output);      /* n_HI/n_H */
  fwrite(temp_H1,sizeof(double),NBINS*NumLos,output);   /* T [K], HI weighted */
  fwrite(veloc_H1,sizeof(double),NBINS*NumLos,output);  /* v_pec [km s^-1], HI weighted */
  fwrite(tau_H1,sizeof(double),NBINS*NumLos,output);    /* HI optical depth */
  fwrite(flux_power,sizeof(float),(NBINS+1)/2*NumLos,output);    /* HI optical depth */
#ifdef HELIUM
  fwrite(n_He2,sizeof(double),NBINS*NumLos,output);     /* n_HeII/n_H */
  fwrite(temp_He2,sizeof(double),NBINS*NumLos,output);   /* T [K], HeII weighted */
  fwrite(veloc_He2,sizeof(double),NBINS*NumLos,output); /* v_pec [km s^-1], HeII weighted */
  fwrite(tau_He2,sizeof(double),NBINS*NumLos,output);   /* HeII optical depth */
#endif
  fwrite(posaxis,sizeof(double),NBINS,output);          /* pixel positions, comoving kpc/h */
  fwrite(velaxis,sizeof(double),NBINS,output);          /* pixel positions, km s^-1 */
  fclose(output);
  FreeLOSMemory();
  return;
}

/*****************************************************************************/
void InitLOSMemory(int NumLos)
{  
  Delta        = (double *) calloc((NumLos * NBINS) , sizeof(double));
  n_H1         = (double *) calloc((NumLos * NBINS) , sizeof(double));
  veloc_H1     = (double *) calloc((NumLos * NBINS) , sizeof(double));
  temp_H1      = (double *) calloc((NumLos * NBINS) , sizeof(double));
  tau_H1       = (double *) calloc((NumLos * NBINS) , sizeof(double));
  posaxis      = (double *) calloc(NBINS , sizeof(double));
  velaxis      = (double *) calloc(NBINS , sizeof(double));
  flux_power   = (float *) calloc(NumLos * (NBINS+1)/2, sizeof(float));
  if(!Delta ||!posaxis  ||   !velaxis || 
   !n_H1 || !veloc_H1 || !temp_H1 || !tau_H1 || !flux_power  )
  {
      fprintf(stderr, "Failed to allocate memory!\n");
      exit(1);
  }
#ifdef HELIUM 
  n_He2         = (double *) calloc((NumLos * NBINS) , sizeof(double)); 
  veloc_He2     = (double *) calloc((NumLos * NBINS) , sizeof(double)); 
  temp_He2      = (double *) calloc((NumLos * NBINS) , sizeof(double)); 
  tau_He2       = (double *) calloc((NumLos * NBINS) , sizeof(double)); 
  if(!n_He2  || !veloc_He2 || ! temp_He2  || ! tau_He2 )
  {
      fprintf(stderr, "Failed to allocate helium memory!\n");
      exit(1);
  }
#endif
}
/*****************************************************************************/

/*****************************************************************************/
void FreeLOSMemory(void)
{  
  free(Delta)     ;
  free(posaxis)   ;
  free(velaxis)   ;
  free(n_H1     ) ;
  free(veloc_H1 ) ;
  free(temp_H1  ) ;
  free(tau_H1   ) ;
  free(flux_power);
#ifdef HELIUM 
  free(n_He2     );
  free(veloc_He2 );
  free(temp_He2  );
  free(tau_He2   );
#endif
}
/*****************************************************************************/
