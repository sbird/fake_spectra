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

/* Physical constants, cgs units */
#define  SIGMA_T 6.652458558e-25 /* Thompson cross-section in cm^2*/
//Note BOLTZMANN is given in cgs units, velocities supplied
//from outside will be here are in km/s
#define  BOLTZMANN  1.3806504e-16  /*  ergs K-1 or cm2 g s-2 K-1 */

#define  LIGHT      2.99792458e10 /*in cm/s*/
// convert energy/unit mass to erg/g
#define  ESCALE 1.0e10
/* Some useful numbers */
#define  GAMMA (5.0/3.0)

/*Conversion factor between internal energy and mu and temperature in K */
#define TSCALE ((GAMMA-1.0) * PROTONMASS * ESCALE / BOLTZMANN)

#define NGRID 8

//Threshold of tau below which we stop computing profiles
#define TAUTAIL 0.001

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
    //return 3./(4*M_PI)*(zhigh - zlow);
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

/* Compute the Voigt or Gaussian profile.
 * Uses the approximation to the Voigt profile from Tepper-Garcia, 2006, MNRAS, 369, 2025
 * (astro-ph/0602124) eq. 25, somewhat simplified
 * which is accurate to 1% in the worst case, a narrow region between the wings and the core.
 * Arguments:
 * T0 = (vdiff/btherm)**2
 * aa: voigt_fac/btherm
 * (note btherm is sqrt(2k T / M))
 */
inline double profile(const double T0, const double aa)
{
    const double T1 = exp(-T0);
  #ifdef VOIGT
    const double T2 = 1.5/T0;
    const double profile_H1 = (T0 < 1.e-6 ? T1 : T1 - aa/sqrt(M_PI)/T0*(T1*T1*(4.0*T0*T0 + 7.0*T0 + 4.0 + T2) - T2 -1.0));
  #else
    const double profile_H1 = T1;
  #endif
    return profile_H1;
}


class SingleAbsorber
{
    public:

        /*This class evaluates the convolution integral for the optical depth.
         * This is:
         *
         * tau = int_{bin} dv' int_{sph kernel support} dv n(v) Phi(v-v')
         * where n is the local density and Phi is the broadening function.
         *
         * Arguments:
         * btherm: thermal b parameter: sqrt(2kT/m)
         * vel: velocity of particle parallel to sightline
         * vdr2: impact parameter from particle center to sightline in velocity units (km/s)
         * vsmooth: smoothing length in velocity units (km/s)
         * aa: voigt_fac/btherm the parameter for the Voight profile
         */
        SingleAbsorber(double bth_i, double vdr2_i, double vsm_i, double aa_i):
            btherm(bth_i), vdr2(vdr2_i), vsmooth(vsm_i), aa(aa_i)
        {};

        /*Evaluate the outer integral in the convolution which gives us the optical depth.
         * This is:
         *
         * tau = int_{bin} dv' tau_kern_inner(v)
         * Arguments:
         * vlow, vhigh: integration limits, which should be v_bin - v_particle
         */
        double tau_kern_outer(const double vlow, const double vhigh)
        {
            double total = tau_kern_inner(vlow)/2.;
            const double deltav=(vhigh-vlow)/NGRID;
            for(int i=1; i<NGRID; ++i)
            {
                const double vv = i*deltav+vlow;
                total += tau_kern_inner(vv);
            }
            total += tau_kern_inner(vhigh)/2.;
            return deltav*total;
        }

    private:
        /*Evaluate the inner integral in the convolution which gives us the optical depth.
         * This is:
         *
         * tau = int_{bin} dv' int_{sph kernel support} dv n(v) Phi(v-v')
         * where n is the local density and Phi is the broadening function.
         *
         * This does the inner function: int{sph kernel support} dv n(v) Phi(v-v')
         * (Strictly speaking internally we use v'' = v - v_part and v''' = v' - v_part)
         * Arguments:
         * vdr2: impact parameter from particle center to sightline in velocity units (km/s)
         * vouter: velocity value for the outer integral. This should be v_bin - v_particle.
         * btherm: thermal b parameter: sqrt(2kT/m)
         * vsmooth: smoothing length in velocity units (km/s)
         */
        inline double tau_kern_inner(const double vouter)
        {
            //Compute SPH kernel support region
            const double vhigh = sqrt(vsmooth*vsmooth-vdr2);
            //Integration region goes from -vhigh to vhigh
            const double deltav=2.*vhigh/NGRID;
            //Because we are integrating over the whole sph kernel,
            //the first and last terms will have q = 1, sph_kernel = 0, so don't need to compute them.
            double total = 0;
            for(int i=1; i<NGRID; ++i)
            {
                const double vv = i*deltav-vhigh;
                const double q = sqrt(vdr2+vv*vv)/vsmooth;
                //The difference between this velocity bin and the particle velocity
                const double vdiff = vv - vouter;
                const double T0 = pow(vdiff/btherm,2);
                total+=sph_kernel(q)*profile(T0,aa);
            }
            return 8*deltav*total/M_PI;
        }

        const double btherm;
        const double vdr2;
        const double vsmooth;
        const double aa;
};



//Factor of 1e5 in bfac converts from cm/s to km/s
LineAbsorption::LineAbsorption(const double lambda, const double gamma, const double fosc, const double amumass, const double velfac_i, const double boxsize, const double atime_i):
sigma_a( sqrt(3.0*M_PI*SIGMA_T/8.0) * lambda  * fosc ),
bfac( sqrt(2.0*BOLTZMANN/(amumass*PROTONMASS))/1e5 ),
voigt_fac( gamma*lambda/(4.*M_PI) ),
velfac(velfac_i), vbox(boxsize*velfac_i), atime(atime_i)
{
}

/* Add the absorption from a particle to the spectrum in the array
 * tau, and the density from the particle to the array colden
 * The slightly C-style interface is so we can easily use the data in python
 */
void LineAbsorption::add_particle(double * tau, double * colden, const int nbins, const double dr2, const float dens, const float ppos, const float pvel, const float temp, const float smooth)
{
  /*Factor to convert the dimensionless quantity found by sph_kern_frac to a column density,
   * in [dens units] * [h units] (atoms/cm^3 * kpc/h comov if from python,
   * (1e10 M_sun /h) / (kpc/h)^2 if from C).
   * The factor of h is because we compute int_z ρ dz, using dimensionless units for z, s.t. χ = z/h,
   */
  const double avgdens = dens * smooth;
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
  // Compute the column density
  for(int z=zlow; z<=zhigh; z++)
  {
      /*Difference between velocity of bin this edge and particle in units of the smoothing length*/
      const double vlow = (boxtosm*z - vel);
      const double vhigh = (boxtosm*(z+1) - vel);

      //colden in units of [den units]*[h units] * integral in terms of z / h
      const double colden_this = avgdens*sph_kern_frac(vlow, vhigh, bb2);

      // The index may be periodic wrapped.
      // Index in units of the box
      int j = z % nbins;
      if (j < 0)
        j+=nbins;
      colden[j] += colden_this;
  }
  //Finish now if not computing absorption
  if (!tau) {
      return;
  }

  /* btherm has the units of velocity: km/s*/
  const double btherm = bfac*sqrt(temp);
  // Create absorption object
  SingleAbsorber absorber ( btherm, velfac*velfac*dr2, vsmooth, voigt_fac/btherm );
  // Do the tau integral for each bin
  const double bintov = vbox/nbins;
  // Amplitude factor for the strength of the transition.
  // sqrt(pi)*c / btherm comes from the profile.
  const double amp = sigma_a / sqrt(M_PI) * (LIGHT/btherm);
  //Bin nearest the particle
  const int zmax = floor(vel/bintov);
  //Go from nearest the particle to half a box away.
  //After that we certainly want to stop as we will have wrapped around.
  //This means that damping wings which span a whole box will be cut off at the edges.
  for(int z=zmax; z<zmax+nbins/2; ++z)
  {
      double vlow = z*bintov - vel;
      const double taulast=amp*absorber.tau_kern_outer(vlow, vlow+bintov);
      tau[z % nbins ]+=taulast;
      //Absorption will only decrease as you go further from the particle.
      if(taulast < TAUTAIL)
        break;
  }
  //Go from the particle backwards
  for(int z=zmax-1; z>zmax-nbins/2; --z)
  {
      double vlow = z*bintov - vel;
      const double taulast=amp*absorber.tau_kern_outer(vlow, vlow+bintov);
      tau[(z + nbins) % nbins ]+=taulast;
      //Absorption will only decrease as you go further from the particle.
      if(taulast < TAUTAIL)
        break;
  }

  return;
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
