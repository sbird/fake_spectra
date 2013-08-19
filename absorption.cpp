/* Copyright (c) 2013 Simeon Bird <spb@ias.edu>
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

#include <cmath>

#define  BOLTZMANN   1.3806504e-23  /* m2 kg s-2 K-1 */
#define  LIGHT           2.99792458e8 /*in km/s*/
#define  PROTONMASS  1.66053886e-27 /* 1 a.m.u */
#define  KPC 3.08568025e19 /*in m*/
#define  SIGMA_T 6.652458558e-29 /* Thompson cross-section in m^2*/

class LineAbsorption
{
    public:
      LineAbsorption(const double lambda, const double gamma, const double fosc, const double amumass):
      sigma_a( sqrt(3.0*M_PI*SIGMA_T/8.0) * lambda  * fosc ),
      bfac( sqrt(2.0*BOLTZMANN/(amumass*PROTONMASS)) ),
      voigt_fac( gamma*lambda/(4.*M_PI) )
      {    }
  
      /* Compute the absorption in a single bin, using 
       * either straight Gaussian or a Voigt profile.
       * Arguments: 
       * colden: column density of absorber in amu per m^2.
       * vdiff: the relative velocities between absorper and bin.
       * temp: temperature of absorber in K
       */
      inline double tau_single(const double colden, const double vdiff, const double temp)
      {
          /* b has the units of velocity: km/s*/
          const double b_H1   = bfac*sqrt(temp);
          const double T0 = pow(vdiff/b_H1,2);
          const double T1 = exp(-T0);
          /* Voigt profile: Tepper-Garcia, 2006, MNRAS, 369, 2025
           * includes thermal and doppler broadening. */
        #ifdef VOIGT
          const double aa_H1 = voigt_fac/b_H1;
          const double T2 = 1.5/T0;
          const double profile_H1 = (T0 < 1.e-6 ? T1 : T1 - aa_H1/sqrt(M_PI)/T0*(T1*T1*(4.0*T0*T0 + 7.0*T0 + 4.0 + T2) - T2 -1.0));
        #else
          const double profile_H1 = T1;
        #endif
          return sigma_a / sqrt(M_PI) * (LIGHT/b_H1) * colden * profile_H1;
      }
  
      /* Absorption cross-sections m^2 */
      const double sigma_a;
      /* Constant factor to turn sqrt(temperature) into velocity*/
      const double bfac;
      /* Factor to turn b into a dimensionless Voigt profile broadening factor, 
       * giving the balance between doppler and thermal broadening. */
      const double voigt_fac;
};

/*****************************************************************************/
/* This function calculates absorption from a given integrated temperature, density
 * and line profile properties.
 * Note: a lot of variables are named _H1. This is historical: the function works for arbitrary species.
 * Arguments are:
 * tau_H1: Array to store the ouput optical depth
 * H1 : species with density, velocity and temperature arrays
 * Hz: conversion factor from linear to velocity space, Hubble(z) in km/s/Mpc
 * box100: box size in comoving kpc/h
 * h100: hubble constant h (~ 0.7)
 * atime: a = 1/(1+z)
 * lambda_lya, gamma_lya, fosc_lya: parameters of the atomic transition (use those from VPFIT)
 * mass: mass of the species in amu
 * */
void Compute_Absorption(double * tau_H1, double * rho, double * veloc, double * temp, const int nbins, const double Hz, const double h100, const double box100, const double atime, const double lambda_lya, const double gamma_lya, const double fosc_lya, const double mass)
{
  /* Conversion factors from internal units */
  const double rscale = (KPC*atime)/h100;   /* convert length to m */
  /*    Calculate the length scales to be used in the box */
  const double vmax = box100 * Hz * rscale/ (1e3*KPC); /* box size (kms^-1) */
  const double dzgrid   = box100 * rscale / (double) nbins; /* bin size m */
  const double dvbin = dzgrid * Hz / (1e3*KPC); /* velocity bin size (kms^-1) */
  LineAbsorption line(lambda_lya, gamma_lya, fosc_lya, mass);
  /* Compute the HI Lya spectra */
  for(int i=0;i<nbins;i++){
      for(int j=0;j<nbins;j++)
        {
          double u_H1, vdiff_H1;

          u_H1  = dvbin*j*1.0e3;
      #ifdef PECVEL
          u_H1 +=veloc[j]*1.0e3;
      #endif
          /* Note this is indexed with i, above with j!
           * This is the difference in velocities between two clouds
           * on the same sightline*/
          vdiff_H1  = fabs(dvbin*i*1.0e3 - u_H1); /* ms^-1 */
       #ifdef PERIODIC
          if (vdiff_H1 > (vmax/2.0*1.0e3))
              vdiff_H1 = (vmax*1.0e3) - vdiff_H1;
       #endif
          tau_H1[j] += line.tau_single(dzgrid * rho[j]/(mass*PROTONMASS), vdiff_H1, temp[j]);
        }
  }             /* Spectrum convolution */
      
  return;
}

