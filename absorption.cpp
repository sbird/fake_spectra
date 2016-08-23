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
#include "singleabs.h"

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

//Threshold of tau below which we stop computing profiles
//Note that because this is for each particle, it should be fairly small.
#define TAUTAIL 1e-5

/* Find the integral of the particle density in this pixel by integrating an SPH kernel
 * over the z direction.
 * Arguments:
 * zlow - Lower z limit for the integral (as z distance from particle center).
 * zhigh - Upper z limit for the integral (again as distance from particle center)
 * smooth - smoothing length
 * dr2 - transverse distance from particle to pixel, squared.
 * zrange - sqrt(smooth*smooth - dr2) (so it can be precomputed)
 *
 * If K(q) is the SPH kernel we need
 * int_zlow^zhigh K(q) dz
 *
 * Normalized such that 4 pi int_{q < 1} K(q) q^2 dq = 1
 * (this matches the definition used in gadget for the smoothing length!)
 * and q^2 = (b^2 + z^2)/h^2, implying that (since int_{q<1} K(q) q^2 dq = 1/32)
 * we want a normalisation of 32/(4*pi)
 *
 * */
double sph_cubic_kern_frac(double zlow, double zhigh, const double smooth, const double dr2, const double zrange)
{
    //Truncate bin size to support of kernel
    zlow = std::max(zlow, -zrange);
    zhigh = std::min(zhigh, zrange);
    if (zlow > zhigh)
        return 0;
    const double qlow = sqrt(dr2+zlow*zlow)/smooth;
    double total = sph_kernel(qlow)/2.;
    const double deltaz=(zhigh-zlow)/NGRID;
    for(int i=1; i<NGRID; ++i)
    {
        const double zz = i*deltaz+zlow;
        const double q = sqrt(dr2+zz*zz)/smooth;

        total+=sph_kernel(q);
    }
    double qhigh = sqrt(dr2+zhigh*zhigh)/smooth;

    total += sph_kernel(qhigh)/2.;
    return deltaz*total;
}

/* Find the fraction of the total particle density in this pixel by integrating a top hat kernel
 * over the z direction. This assumes that rho = rho_0, a constant, within the cell.
 * Arguments:
 * zlow - Lower z limit for the integral (as z distance from particle center).
 * zhigh - Upper z limit for the integral (again as distance from particle center)
 * smooth - smoothing length
 * dr2 - transverse distance from particle to pixel, squared.
 * zrange - sqrt(smooth*smooth - dr2) (so it can be precomputed)
 * */

double tophat_kern_frac(double zlow, double zhigh, const double smooth, const double dr2, const double zrange)
{
    //Integration limits
    zlow = std::max(zlow, -zrange);
    zhigh = std::min(zhigh, zrange);
    return std::max(0.,zhigh - zlow);
}

kern_frac_func get_kern_frac(const int kernel)
{
    if(kernel == TOP_HAT_KERNEL)
        return tophat_kern_frac;
    else
        return sph_cubic_kern_frac;
}

//Factor of 1e5 in bfac converts from cm/s to km/s
//Factor of 1e5 in voigt_fac converts from cm/s to km/s
LineAbsorption::LineAbsorption(const double lambda, const double gamma, const double fosc, const double amumass, const double velfac_i, const double boxsize, const double atime_i,const int kernel_i):
sigma_a( sqrt(3.0*M_PI*SIGMA_T/8.0) * lambda  * fosc ),
bfac( sqrt(2.0*BOLTZMANN/(amumass*PROTONMASS))/1e5 ),
voigt_fac( gamma*lambda/(4.*M_PI)/1e5 ),
velfac(velfac_i), vbox(boxsize*velfac_i), atime(atime_i),
kern_frac(get_kern_frac(kernel_i))
{
}

/* Add the absorption from a particle to the spectrum in the array
 * tau, and the density from the particle to the array colden
 * The slightly C-style interface is so we can easily use the data in python
 */
void LineAbsorption::add_colden_particle(double * colden, const int nbins, const double dr2, const float dens, const float pos, const float smooth)
{
  //If we are outside the kernel, do nothing.
  if (smooth*smooth - dr2 <= 0)
      return;
  //z range covered by particle in kpc/h
  const double zrange = sqrt(smooth*smooth - dr2);
  //Conversion between units of position to units of the box.
  const double boxtokpc = vbox / nbins / velfac;
  // z is position in units of the box
  const int zlow = floor((pos - zrange) / boxtokpc);
  const int zhigh = ceil((pos + zrange) / boxtokpc);
  // Compute the column density
  for(int z=zlow; z<=zhigh; z++)
  {
      //Difference between position of bin this edge and particle
      const double plow = (boxtokpc*z - pos);
      // The index may be periodic wrapped.
      // Index in units of the box
      int j = z % nbins;
      if (j<0)
        j+=nbins;
      /*Factor to convert the quantity found by sph_kern_frac which has units of kpc/h, to a column density,
       * in [dens units] * [h units] (atoms/cm^3 * kpc/h if from python,
       * (1e10 M_sun /h) / (kpc/h)^2 if from C).
       * We compute int_z ρ dz
       */
      //colden in units of [den units]*[h units] * integral in terms of z
      colden[j] += dens*kern_frac(plow, plow + boxtokpc, smooth, dr2,zrange);
  }
}

void LineAbsorption::add_tau_particle(double * tau, const int nbins, const double dr2, const float dens, const float ppos, const float pvel, const float temp, const float smooth)
{
  /* Velocity of particle parallel to los: pos in kpc/h comoving
     to vel in km/s physical. Note that gadget velocities come comoving,
     so we need the sqrt(a) conversion factor.
   */
  const double vel = velfac * ppos + pvel * sqrt(atime);
  /* btherm has the units of velocity: km/s*/
  const double btherm = bfac*sqrt(temp);
  //Double check we are within the kernel support
  if(smooth*smooth - dr2 <=0)
      return;
  // Create absorption object
  SingleAbsorber absorber ( btherm, velfac*velfac*dr2, velfac*smooth, voigt_fac/btherm );
  // Do the tau integral for each bin
  const double bintov = vbox/nbins;
  // Amplitude factor for the strength of the transition.
  // sqrt(pi)*c / btherm comes from the profile.
  //Extra factor of 1e5 because LIGHT is in cm/s and btherm is in km/s
  const double amp = sigma_a / sqrt(M_PI) * (LIGHT/1e5/btherm);
  //Bin nearest the particle
  const int zmax = floor(vel/bintov);
  //Go from nearest the particle to half a box away.
  //After that we certainly want to stop as we will have wrapped around.
  //This means that damping wings which span a whole box will be cut off at the edges.
  for(int z=zmax; z<zmax+nbins/2; ++z)
  {
      double vlow = z*bintov - vel;
      //dens is in atoms/cm^3, amp is in cm^2. tau_kern_outer returns in velocity units (here km/s),
      //so then we have units of the end of km/s /cm, and need to divide by velfac to get kpc/h/cm
      const double taulast=amp*dens*absorber.tau_kern_outer(vlow, vlow+bintov)/velfac;
      // Make sure index is properly wrapped
      int j = z % nbins;
      if (j<0)
        j+=nbins;
      tau[j]+=taulast;
      //Absorption will only decrease as you go further from the particle.
      if(taulast < TAUTAIL)
        break;
  }
  //Go from the particle backwards
  for(int z=zmax-1; z>=zmax-nbins/2; --z)
  {
      double vlow = z*bintov - vel;
      const double taulast=amp*dens*absorber.tau_kern_outer(vlow, vlow+bintov)/velfac;
      // Make sure index is properly wrapped
      int j = z % nbins;
      if (j<0)
        j+=nbins;
      tau[j]+=taulast;
      //Absorption will only decrease as you go further from the particle.
      if(taulast < TAUTAIL)
        break;
  }
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
