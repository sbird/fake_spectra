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
//Note that because this is for each particle, it should be fairly small.
#define TAUTAIL 1e-5

inline double sph_kernel(const double q)
{
    if (q >= 1)
        return 0;
    if(q<0.5)
        return 1-6*q*q+6*q*q*q;
    else
        return 2*pow(1.- q,3);
}

#ifndef TOP_HAT_KERNEL
/* Find the integral of the particle density in this pixel by integrating an SPH kernel
 * over the z direction.
 * Arguments:
 * zlow - Lower z limit for the integral (as z distance from particle center).
 * zhigh - Upper z limit for the integral (again as distance from particle center)
 * bb2 - transverse distance from particle to pixel, squared.
 *
 * If K(q) is the SPH kernel we need
 * int_zlow^zhigh K(q) dz
 *
 * Normalized such that 4 pi int_{q < 1} K(q) q^2 dq = 4 pi/3
 * and q^2 = (b^2 + z^2)/h^2, implying that (since int_{q<1} K(q) q^2 dq = 1/32)
 * we want a normalisation of 32/3
 *
 * */
double sph_kern_frac(double zlow, double zhigh, double smooth, double dr2)
{
    //Outside useful range.
    const double zrange = sqrt(smooth*smooth - dr2);
    if (zlow > zrange || zhigh < -zrange){
        return 0;
    }
    //Maximal range that will do anything
    zlow = std::max(zlow, -zrange);
    zhigh = std::min(zhigh, zrange);
    double total = sph_kernel(sqrt(dr2+zlow*zlow)/smooth)/2.;
    const double deltaz=(zhigh-zlow)/NGRID;
    for(int i=1; i<NGRID; ++i)
    {
        const double zz = i*deltaz+zlow;
        const double q = sqrt(dr2+zz*zz)/smooth;
        total+=sph_kernel(q);
    }
    double qhigh = sqrt(dr2+zhigh*zhigh)/smooth;
    total += sph_kernel(qhigh)/2.;
    return 32./3.*deltaz*total;
}

#else

/* Find the fraction of the total particle density in this pixel by integrating a top hat kernel
 * over the z direction. This assumes that rho = rho_0, a constant, within the cell.
 * Arguments:
 * zlow - Lower z limit for the integral (as z distance from particle center).
 * zhigh - Upper z limit for the integral (again as distance from particle center)
 * smooth - smoothing length
 * dr2 - transverse distance from particle to pixel, squared.
 * */

double sph_kern_frac(double zlow, double zhigh, double smooth, double dr2)
{
    //Cell boundaries
    const double zrange = sqrt(smooth*smooth - dr2);
    //Integration limits
    zlow = std::max(zlow, -zrange);
    zhigh = std::min(zhigh, zrange);
    return std::max(0.,zhigh - zlow);
}
#endif

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
    const double profile_H1 = (T0 < 1.e-6 ? T1 : T1 - aa/sqrt(M_PI)/T0*(T1*T1*(4.0*T0*T0 + 7.0*T0 + 4.0 + 1.5/T0) - 1.5/T0 -1.0));
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

        /*Find the mean optical depth in the bin, as bins may be broader than the SPH kernel.
         * This is:
         *
         * tau = 1/V int_{bin} dv' tau_kern_inner(v)
         * where V is the bin width
         * Arguments:
         * vlow, vhigh: integration limits, which should be v_bin - v_particle
         */
        double tau_kern_outer(const double vlow, const double vhigh)
        {
            double total = tau_kern_inner(vlow)/2.;
            //Note that the number of points does not need to be NGRID here,
            //but it is for now.
            const double deltav=(vhigh-vlow)/NGRID;
            for(int i=1; i<NGRID; ++i)
            {
                const double vv = i*deltav+vlow;
                total += tau_kern_inner(vv);
            }
            total += tau_kern_inner(vhigh)/2.;
            return total/NGRID;
        }

    private:
        /* This gives us the optical depth at velocity vouter, from a convolution of the density and
         * broadening function.
         * This is:
         *
         * tau(v') = int_{sph kernel support} dv n(v) Phi(v-v')
         * where n is the local density and Phi is the broadening function.
         *
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
            return 32./3.*deltav*total;
        }

        const double btherm;
        const double vdr2;
        const double vsmooth;
        const double aa;
};



//Factor of 1e5 in bfac converts from cm/s to km/s
//Factor of 1e5 in voigt_fac converts from cm/s to km/s
LineAbsorption::LineAbsorption(const double lambda, const double gamma, const double fosc, const double amumass, const double velfac_i, const double boxsize, const double atime_i):
sigma_a( sqrt(3.0*M_PI*SIGMA_T/8.0) * lambda  * fosc ),
bfac( sqrt(2.0*BOLTZMANN/(amumass*PROTONMASS))/1e5 ),
voigt_fac( gamma*lambda/(4.*M_PI)/1e5 ),
velfac(velfac_i), vbox(boxsize*velfac_i), atime(atime_i)
{
}

/* Add the absorption from a particle to the spectrum in the array
 * tau, and the density from the particle to the array colden
 * The slightly C-style interface is so we can easily use the data in python
 */
void LineAbsorption::add_colden_particle(double * colden, const int nbins, const double dr2, const float dens, const float ppos, const float pvel, const float smooth)
{
  /* Velocity of particle parallel to los: pos in kpc/h comoving
     to vel in km/s physical. Note that gadget velocities come comoving,
     so we need the sqrt(a) conversion factor.*/
  const double pos = ppos + pvel * sqrt(atime)/velfac;
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
      colden[j] += dens*sph_kern_frac(plow, plow + boxtokpc, smooth, dr2);
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
  for(int z=zmax-1; z>zmax-nbins/2; --z)
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
