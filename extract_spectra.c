
#include "headers.h"
#include "global_vars.h"
#include "parameters.h"

double *rhoker_H,*Delta,*posaxis,*velaxis;
double *rhoker_H1,*velker_H1,*temker_H1;
double *rhoker_He2,*velker_He2,*temker_He2;
double *n_H1,*veloc_H1,*temp_H1,*tau_H1;
double *n_He2,*veloc_He2,*temp_He2,*tau_He2;


/*****************************************************************************/
void SPH_interpolation()
{
  double Hz,rhoc,critH,H0,mu,rscale,vscale,mscale,escale,hscale;
  double zmingrid,zmaxgrid,dzgrid,dzinv,boxsize,box2,dzbin,vmax,dvbin;
  double xproj,yproj,zproj,xx,yy,zz,hh,h2,h4,dr,dr2;
  double hinv2,hinv3,vr,Temperature,dzmax,zgrid;
  double deltaz,dz,dist2,H1frac,He2frac,q,kernel,velker,temker;
  double sigma_Lya_H1,sigma_Lya_He2,vmax2,vdiff_H1,vdiff_He2;
  double A_H1,T0,T1,T2,tau_H1j,aa_H1,u_H1,b_H1,profile_H1;
  double A_He2,T3,T4,T5,tau_He2j,aa_He2,u_He2,b_He2,profile_He2;
  int i,iproc,ic,iaxis,iz,ioff,j,iiz,ii,jj;
  
  double ztime[1];

  FILE *output;
  

  InitLOSMemory();
  printf("Allocating memory...done\n");
  
  srand48(241008); /* random seed generator */
  
  H0 = 1.0e5/MPC; /* 100kms^-1Mpc^-1 in SI */ 
  Hz = 100.0*h100 * sqrt(1.+omega0*(1./atime-1.)+omegaL*((atime*atime) -1.))/atime;  
  
  /* Critical matter/energy density at z = 0.0 */
  rhoc = 3.0 * (H0*h100)*(H0*h100) / (8.0 * PI * GRAVITY); /* kgm^-3 */
  
  /* Mean hydrogen mass density of the Universe */
  critH = (rhoc * OMEGAB * XH) / (atime*atime*atime); /* kgm^-3*/
  
  /* Absorption cross-sections m^2 */
  sigma_Lya_H1  = sqrt(3.0*PI*SIGMA_T/8.0) * LAMBDA_LYA_H1  * FOSC_LYA;
  /* sigma_Lya_He2 = sqrt(3.0*PI*SIGMA_T/8.0) * LAMBDA_LYA_HE2 * FOSC_LYA; */
  
  /* Conversion factors from internal units */
  rscale = (KPC*atime)/h100;   /* convert length to m */
  vscale = sqrt(atime);        /* convert velocity to kms^-1 */
  mscale = (1.0e10*SOLAR_MASS)/h100; /* convert mass to kg */
  escale = 1.0e6;           /* convert energy/unit mass to J kg^-1 */
  hscale = rscale * 0.5; /* Note the factor of 0.5 for this kernel definition */
  
  /*   Convert to SI units from GADGET-3 units */
  for(i=1;i<NumPart[0]+1;i++)
    {
      for(ic=0;ic<3;ic++)
	{
	  P[i].Pos[ic] *= rscale; /* m, physical */
	  P[i].Vel[ic] *= vscale; /* km s^-1, physical */
	}
      
      P[i].h *= hscale;   /* m, physical */
      P[i].Mass_d = (double)P[i].Mass * mscale;   /* kg */
      //   printf("try to print mass %e %e \n",mscale,P[i].Mass_d);

      /* Mean molecular weight */
      mu = 1.0/(XH*(0.75+P[i].Ne) + 0.25);
      P[i].U *= ((GAMMA-1.0) * mu * HMASS * PROTONMASS * escale ) / BOLTZMANN; /* K */
    }
  printf("Converting units...done\n");
  
  
  /*    Calculate the length scales to be used in the box */
  zmingrid = 0.0;
  zmaxgrid = box100*rscale;  /* box sizes in physical m */
  dzgrid   = (zmaxgrid-zmingrid) / (double)NBINS; /* bin size (physical m) */
  dzinv    = 1. / dzgrid;
  boxsize  = zmaxgrid;   
  box2     = 0.5 * boxsize;
  dzbin = box100/ (double)NBINS; /* bin size (comoving kpc/h) */
  
  vmax = box100 * Hz * rscale/ MPC; /* box size (kms^-1) */
  vmax2 = vmax/2.0; /* kms^-1 */
  dvbin = vmax / (double)NBINS; /* bin size (kms^-1) */
  
 

  /*   Initialise distance coordinate for iaxis */
  posaxis[0]=0.0;
  velaxis[0]=0.0;
  
  for(i=0;i<NBINS-1;i++)
    {
      posaxis[i+1] = posaxis[i] + dzbin; /* comoving kpc/h */
      velaxis[i+1] = velaxis[i] + dvbin; /* physical km s^-1 */
    }
  
  
  /*    Generate random coordinates for a point in the box */
  for(iproc=0;iproc<NUMLOS;iproc++)
    { 
      /*Pick a random sightline*/
      do	
      	iaxis = (int)(drand48()*4);
      while (iaxis == 0 || iaxis==4); 
      
      xproj = drand48()*box100*rscale;
      yproj = drand48()*box100*rscale;
      zproj = drand48()*box100*rscale;
      
      printf("Interpolating line of sight %d...done\n",iproc);
      
      /* Loop over particles in LOS and do the SPH interpolation */
      /* This first finds which particles are near this sight line. 
       * Probably a faster way to do that. 
       * Then adds the total density, temp. and velocity for near particles to 
       * the binned totals for that sightline*/
      for(i=1;i<NumPart[0]+1;i++)
	{
	  
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
          else if (iaxis == 3)
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
		   /* He2frac = P[i].NHep; *//* nHeII/nH */
		   
		   /* Central vertex to contribute to */
		   if (iaxis == 1)
		     iz = (xx - zmingrid) * dzinv +1  ;
		   else if (iaxis == 2) 
		     iz = (yy - zmingrid) * dzinv +1 ;
		   else if (iaxis == 3)
		     iz = (zz - zmingrid) * dzinv +1;
		   
		   dzmax = sqrt(fabs(h4 - dr2));
		   ioff = (int)(dzmax * dzinv) +1;
		   
		   /* Loop over contributing vertices */
		   for(iiz = iz-ioff; iiz < iz+ioff+1 ; iiz++)
		     {
		       j = iiz;
		       j = ((j-1+10*NBINS) % NBINS);
		       
		       zgrid = zmingrid + (double)(j) * dzgrid;
		       
		      if (iaxis == 1)
                        deltaz = zgrid - xx;
		      else if (iaxis == 2)
                        deltaz = zgrid - yy;
		      else if (iaxis == 3)
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
			    kernel = (1.+ (q*q) * (-1.5 + 0.75 * q) )/PI;
			  else
			    kernel = 0.25*(2.0-q)*(2.0-q)*(2.0-q)/PI;
			  
			  kernel = kernel * hinv3;  

			  kernel = P[i].Mass_d * kernel; /* kg m^-3 */
			  velker = vr * kernel; /* kg m^-3 * km s^-1 */
			  temker = Temperature * kernel; /* kg m^-3 * K */
			  
			  jj = j + ((NBINS)*iproc);
			  rhoker_H[jj]  += kernel * XH;		 
			  
			  rhoker_H1[jj] += kernel * XH * H1frac;
			  velker_H1[jj] += velker * XH * H1frac;
			  temker_H1[jj] += temker * XH * H1frac;

			  //printf("try %d %e.. %e .done\n",jj,kernel,P[i].Mass_d);		      

			  
			}      /* dist2 < 4h^2 */
		     }        /* loop over contributing vertices */
		 }           /* dx^2+dy^2 < 4h^2 */
	    }               /* dx < 2h */
	}                  /* Loop over particles in LOS */
      
      
      for(i = 0;i<NBINS;i++)
	{
	  ii = i + (NBINS*iproc);
	  
	  Delta[ii]     = log10(rhoker_H[ii]/critH);   /* log H density normalised by mean 
                                                          H density of universe */
	  n_H1[ii]      = rhoker_H1[ii]/rhoker_H[ii];  /* HI/H */
	  veloc_H1[ii]  = velker_H1[ii]/rhoker_H1[ii]; /* HI weighted km s^-1 */ 
	  temp_H1[ii]   = temker_H1[ii]/rhoker_H1[ii]; /* HI weighted K */
      	}
      
      /* Prefactor for optical depth  */
      A_H1 = sigma_Lya_H1*C*dzgrid/sqrt(PI);
      A_He2 = sigma_Lya_He2*C*dzgrid/sqrt(PI);
    
      
      /* Compute the HI and He2 Lya spectra */
      for(i=0;i<NBINS;i++)
	{
	  for(j=0;j<NBINS;j++)
	    {
	      
	      jj =  j + (NBINS*iproc);
	      
	      if (PECVEL == 1)
		{
		  u_H1 =  (velaxis[j] + veloc_H1[jj])*1.0e3;
		  /* u_He2 = (velaxis[j] + veloc_He2[jj])*1.0e3; */
		}
	      else
		{
		  u_H1  = velaxis[j]*1.0e3;
		  /*  u_He2 = velaxis[j]*1.0e3; */
		}
		  
              /* Note this is indexed with i, above with j! 
               * This is the difference in velocities between two clouds 
               * on the same sightline*/
	      vdiff_H1  = fabs((velaxis[i]*1.0e3) - u_H1); /* ms^-1 */
	      /* vdiff_He2 = fabs((velaxis[i]*1.0e3) - u_He2);*/ /* ms^-1 */
	      
	      
	      if (PERIODIC == 1)
		{
		  if (vdiff_H1 > (vmax2*1.0e3))
		    vdiff_H1 = (vmax*1.0e3) - vdiff_H1;
		  
		  /* if (vdiff_He2 > (vmax2*1.0e3))
		     vdiff_He2 = (vmax*1.0e3) - vdiff_He2; */
		}
	
	      
	      b_H1   = sqrt(2.0*BOLTZMANN*temp_H1[jj]/(HMASS*PROTONMASS));
	      b_He2  = sqrt(2.0*BOLTZMANN*temp_He2[jj]/(HEMASS*PROTONMASS));
	      
	      T0 = (vdiff_H1/b_H1)*(vdiff_H1/b_H1);
	      T1 = exp(-T0);
	      
	      T3 = (vdiff_He2/b_He2)*(vdiff_He2/b_He2);
	      T4 = exp(-T3);
	      
	      /* Voigt profile: Tepper-Garcia, 2006, MNRAS, 369, 2025 */ 
	      if (VOIGT == 1)
		{
		  aa_H1 = GAMMA_LYA_H1*LAMBDA_LYA_H1/(4.0*PI*b_H1);
		  aa_He2 = GAMMA_LYA_HE2*LAMBDA_LYA_HE2/(4.0*PI*b_He2);
		  
		  T2 = 1.5/T0;	
		  T5 = 1.5/T3; 
		  
		  if(T0 < 1.0e-6)
		    profile_H1  = T1;
		  else
		    profile_H1  = T1 - aa_H1/sqrt(PI)/T0 
		      *(T1*T1*(4.0*T0*T0 + 7.0*T0 + 4.0 + T2) - T2 -1.0);
		  
		  /* if(T3 < 1.0e-6)
		       profile_He2  = T4;
 		  else
		    profile_He2 = T4 - aa_He2/sqrt(PI)/T3 
		    *(T4*T4*(4.0*T3*T3 + 7.0*T3 + 4.0 + T5) - T5 -1.0); */
		  
		}
	      
	      else
		{
		  profile_H1 = T1;
		  /* profile_He2 = T4; */
		}
	      
	      
	      tau_H1j  = A_H1  * rhoker_H1[jj]  * profile_H1  /(HMASS*PROTONMASS*b_H1);
	      /*   tau_He2j = A_He2 * rhoker_He2[jj] * profile_He2 /(HMASS*PROTONMASS*b_He2); */
	      
	      ii =  i + (NBINS*iproc);
	      
	      tau_H1[ii]  += tau_H1j;
	      /* tau_He2[ii] += tau_He2j; */
	      // printf("pixel %d.. %lf .done\n",ii,T1);
	    }
	}             /* Spectrum convolution */
    }                /* Loop over numlos random LOS */
  
  ztime[0] = 1.0/atime - 1.0;
  
  output = fopen("spectra1024.dat","wb");
  fwrite(ztime,sizeof(double),1,output);
  fwrite(Delta,sizeof(double),NBINS*NUMLOS,output);     /* gas overdensity */
  fwrite(n_H1,sizeof(double),NBINS*NUMLOS,output);      /* n_HI/n_H */
  fwrite(temp_H1,sizeof(double),NBINS*NUMLOS,output);   /* T [K], HI weighted */
  fwrite(veloc_H1,sizeof(double),NBINS*NUMLOS,output);  /* v_pec [km s^-1], HI weighted */
  fwrite(tau_H1,sizeof(double),NBINS*NUMLOS,output);    /* HI optical depth */
  /*  fwrite(n_He2,sizeof(double),NBINS*NUMLOS,output); */    /* n_HeII/n_H */
  /* fwrite(temp_He2,sizeof(double),NBINS*NUMLOS,output); */  /* T [K], HeII weighted */
  /* fwrite(veloc_He2,sizeof(double),NBINS*NUMLOS,output);*/ /* v_pec [km s^-1], HeII weighted */
  /* fwrite(tau_He2,sizeof(double),NBINS*NUMLOS,output);*/   /* HeII optical depth */
  fwrite(posaxis,sizeof(double),NBINS,output);          /* pixel positions, comoving kpc/h */
  fwrite(velaxis,sizeof(double),NBINS,output);          /* pixel positions, km s^-1 */
  fclose(output);
  
  return;
}

/*****************************************************************************/
void InitLOSMemory(void)
{  
  rhoker_H     = (double *) malloc((NUMLOS * NBINS) * sizeof(double));
  Delta        = (double *) malloc((NUMLOS * NBINS) * sizeof(double));
  
  posaxis      = (double *) malloc(NBINS * sizeof(double));
  velaxis      = (double *) malloc(NBINS * sizeof(double));
  
  rhoker_H1    = (double *) malloc((NUMLOS * NBINS) * sizeof(double));
  velker_H1    = (double *) malloc((NUMLOS * NBINS) * sizeof(double));
  temker_H1    = (double *) malloc((NUMLOS * NBINS) * sizeof(double));
  
  n_H1         = (double *) malloc((NUMLOS * NBINS) * sizeof(double));
  veloc_H1     = (double *) malloc((NUMLOS * NBINS) * sizeof(double));
  temp_H1      = (double *) malloc((NUMLOS * NBINS) * sizeof(double));
  tau_H1       = (double *) malloc((NUMLOS * NBINS) * sizeof(double));
  
  /* rhoker_He2    = (double *) malloc((NUMLOS * NBINS) * sizeof(double)); */
  /* velker_He2    = (double *) malloc((NUMLOS * NBINS) * sizeof(double)); */
  /* temker_He2    = (double *) malloc((NUMLOS * NBINS) * sizeof(double)); */
  
  /* n_He2         = (double *) malloc((NUMLOS * NBINS) * sizeof(double)); */
  /* veloc_He2     = (double *) malloc((NUMLOS * NBINS) * sizeof(double)); */
  /* temp_He2      = (double *) malloc((NUMLOS * NBINS) * sizeof(double)); */
  /* tau_He2       = (double *) malloc((NUMLOS * NBINS) * sizeof(double)); */
}
/*****************************************************************************/
