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

#include "absorption.h"
#include <algorithm>

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

#define NGRID 8

inline double sph_kernel(const double q)
{
    if (q >= 1)
        return 0;
    if(q<0.5)
        return 1-6*q*q+6*q*q*q;
    else
        return 2*pow(1.- q,3);
}

/* Find the fraction of the total particle density in this pixel by integrating an SPH kernel
 * over the z direction. All distances in units of smoothing length.
 * Arguments:
 * zlow - Lower z limit for the integral (as z distance from particle center).
 * zhigh - Upper z limit for the integral (again as distance from particle center)
 * bb2 - transverse distance from particle to pixel, squared.
 *
 * If K(q) is the SPH kernel normalized such that 4 pi int_{q < 1} K(q) q^2 dq = 1
 * and q^2 = b^2 + z^2, then this is:
 * int_zlow^zhigh K(q) dz
 * */
double sph_kern_frac(double zlow, double zhigh, double bb2)
{
    //Outside useful range.
    if (zlow > sqrt(1-bb2) || zhigh < -sqrt(1-bb2)){
        return 0;
    }
    //Maximal range that will do anything
    zlow = std::max(zlow, -sqrt(1-bb2));
    zhigh = std::min(zhigh, sqrt(1-bb2));
    double total = sph_kernel(sqrt(bb2+zlow*zlow))/2.;
    const double deltaz=(zhigh-zlow)/NGRID;
    for(int i=1; i<NGRID; ++i)
    {
        const double zz = i*deltaz+zlow;
        const double q = sqrt(bb2+zz*zz);
        total+=sph_kernel(q);
    }
    double qhigh = sqrt(bb2+zhigh*zhigh);
    total += sph_kernel(qhigh)/2.;
    return 8*deltaz*total/M_PI;
}


/* Add the absorption from a particle to the spectrum in the array
 * tau, and the density from the particle to the array colden
 * The slightly C-style interface is so we can easily use the data in python
 */
void ComputeLineAbsorption::add_particle(double * tau, double * colden, const int nbins, const double dr2, const float mass, const float ppos, const float pvel, const float temp, const float smooth)
{
  /*Factor to convert the dimensionless quantity found by sph_kern_frac to a column density,
   * in (1e10 M_sun /h) / (kpc/h)^2.
   * We compute int_z ρ dz, so the correct dimensional factors are:
   * δz M/h^3: bin width * density */
  const double binwdth = (vbox/velfac/nbins);
  const double avgdens = binwdth*mass*pow(smooth,-3);
  /*Impact parameter in units of the smoothing length */
  const double bb2 = dr2/smooth/smooth;
  /* Velocity of particle parallel to los: pos in */
  const double vel = (velfac * ppos + pvel );
  const double vsmooth = velfac * smooth;
  const double zrange = velfac * sqrt(smooth*smooth - dr2);
  const int zlow = floor((nbins/vbox) * (vel - zrange));
  const int zhigh = ceil((nbins/vbox) * (vel + zrange));
  /* Compute the HI Lya spectra */
  for(int z=zlow; z<=zhigh; z++)
  {
      /*Difference between velocity of bin this edge and particle*/
      const double vlow = (vbox*z/nbins - vel)/vsmooth;
      const double vhigh = (vbox*(z+1)/nbins - vel)/vsmooth;

      const double colden_this = avgdens*sph_kern_frac(vlow, vhigh, bb2);

      // The index may be periodic wrapped.
      int j = z % nbins;
      if (j < 0)
        j+=nbins;
      colden[j] += colden_this;
      /* Loop again, because the column density that contributes to this
       * bin may be broadened thermal or doppler broadened*/
      //Add natural broadening someday
      // 
      if (tau) {
        for(int i=0;i<nbins;i++)
        {
            double vdiff = fabs(vbox*(i-j)/nbins);
            if (vdiff > (vbox/2.0))
              vdiff = vbox - vdiff;
            tau[i] += tau_single(colden_this, vdiff, temp);
        }
      }
  }

  return;
}
