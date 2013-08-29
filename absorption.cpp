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
#include <cmath>

/* Physical constants, SI units */
#define  SIGMA_T 6.652458558e-29 /* Thompson cross-section in m^2*/
#define  BOLTZMANN   1.3806504e-23  /* m2 kg s-2 K-1 */
#define  LIGHT           2.99792458e8 /*in m/s*/
// convert energy/unit mass to J kg^-1
#define  ESCALE (1.0e6)
/* Some useful numbers */
#define  GAMMA (5.0/3.0)

/*Conversion factor between internal energy and mu and temperature in K */
#define TSCALE ((GAMMA-1.0) * PROTONMASS * ESCALE / BOLTZMANN)

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

LineAbsorption::LineAbsorption(const double lambda, const double gamma, const double fosc, const double amumass, const double velfac_i, const double boxsize, const double atime_i):
sigma_a( sqrt(3.0*M_PI*SIGMA_T/8.0) * lambda  * fosc ),
bfac( sqrt(2.0*BOLTZMANN/(amumass*PROTONMASS)) ),
voigt_fac( gamma*lambda/(4.*M_PI) ),
velfac(velfac_i), vbox(boxsize*velfac_i), atime(atime_i)
{
}

/* Add the absorption from a particle to the spectrum in the array
 * tau, and the density from the particle to the array colden
 * The slightly C-style interface is so we can easily use the data in python
 */
void LineAbsorption::add_particle(double * tau, double * colden, const int nbins, const double dr2, const float mass, const float ppos, const float pvel, const float temp, const float smooth)
{
  /*Factor to convert the dimensionless quantity found by sph_kern_frac to a column density,
   * in (1e10 M_sun /h) / (kpc/h)^2.
   * We compute int_z ρ dz, using dimensionless units for z, s.t. χ = z/h,
   * ρ = M/V  and h^3 = 3 V /(4 π)
   * so the correct dimensional factors are:
   *  3/(4π) M/h^2 */
  const double avgdens = 3/4./M_PI*mass*pow(smooth,-2);
  /*Impact parameter in units of the smoothing length */
  const double bb2 = dr2/smooth/smooth;
  const double vsmooth = velfac * smooth;
  /* Velocity of particle parallel to los: pos in kpc/h comoving
     to vel in km/s physical. Note that gadget velocities come comoving,
     so we need the sqrt(a) conversion factor.
     Finally divide by h * velfac to give the velocity in units of the smoothing length.*/
  const double vel = (velfac * ppos + pvel * sqrt(atime))/vsmooth;
  //Allowed z range in units of smoothing length
  const double zrange = sqrt(1. - bb2);
  //Conversion between units of the smoothing length to units of the box.
  const double boxtosm = vbox / vsmooth / nbins;
  // z is position in units of the box
  const int zlow = floor((vel - zrange) / boxtosm);
  const int zhigh = ceil((vel + zrange) / boxtosm);
  /* Compute the HI Lya spectra */
  for(int z=zlow; z<=zhigh; z++)
  {
      /*Difference between velocity of bin this edge and particle in units of the smoothing length*/
      const double vlow = (boxtosm*z - vel);
      const double vhigh = (boxtosm*(z+1) - vel);

      //colden in units of Gadget_mass / Gadget_length^2 * integral in terms of z / h
      const double colden_this = avgdens*sph_kern_frac(vlow, vhigh, bb2);

      // The index may be periodic wrapped.
      // Index in units of the box
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

inline double LineAbsorption::tau_single(const double colden, const double vdiff, const double temp)
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

/* Compute temperature (in K) from internal energy.
 * uu: internal energy in Gadget units
 * ne: electron abundance
 * xh: hydrogen mass fraction (0.76)
 * Factor to convert U (J/kg) to T (K) : U = N k T / (γ - 1)
 * T = U (γ-1) μ m_P / k_B
 * where k_B is the Boltzmann constant
 * γ is 5/3, the perfect gas constant
 * m_P is the proton mass
 * μ is 1 / (mean no. molecules per unit atomic weight) calculated in loop.
 */
double compute_temp(const double uu, const double ne, const double xh)
{
    /*Mean molecular weight:
     * \mu = 1 / molecules per unit atomic weight
     *     = 1 / (X + Y /4 + E)
     *     where E = Ne * X, and Y = (1-X).
     *     Can neglect metals as they are heavy.
     *     Leading contribution is from electrons, which is already included
     *     [+ Z / (12->16)] from metal species
     *     [+ Z/16*4 ] for OIV from electrons. */
    const double mu = 1.0/(xh*(0.75+ne) + 0.25);
    return uu*mu * TSCALE; /* T in K */

}
