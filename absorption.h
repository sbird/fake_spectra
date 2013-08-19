#ifndef ABSORPTION_H
#define ABSORPTION_H

#include <cmath>
#include "parameters.h"

/* Class to compute the absorption in a spectral bin from particles. Base class, extended by ComputeLineAbsorption.*/
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

/* Class which extends LineAbsorption to allow computation of absorption from a single particle
 * onto a whole spectrum, instead of just one bin
 * */
class ComputeLineAbsorption: public LineAbsorption
{
    public:
        /*Dimensions:
         * velocities in km/s (physical).
         * distances in kpc/h (comoving)
         * velfac: factor to convert from distance to velocity units.
         * Should be h100 * atime * Hz/1e3 (Hz in km/s/Mpc)
         * */
        //This makes dvbin be in km/s: the 1e3 converts Hz from km/s/Mpc to km/s/kpc
        // kpc /h       h                    km/s/kpc
        //vbox = ( box100 * h100 * atime * Hz/1e3 ) /* box size (kms^-1) */
        ComputeLineAbsorption(const double lambda, const double gamma, const double fosc, const double mass, const double velfac_i, const double boxsize):
        LineAbsorption(lambda, gamma, fosc, mass),
        amumass(mass), velfac(velfac_i), vbox(boxsize*velfac_i)
        {
        }

        /* Add the absorption from a particle to the spectrum in the array
         * tau, and the density from the particle to the array colden
         * The slightly C-style interface is so we can easily use the data in python.
         *
         * Output:
         * tau: array specifying the optical depth of the spectrum.
         * If this is NULL, just compute the column density.
         * colden: array specifying the column density of the spectrum. (atoms/ (comoving kpc/h)^2)
         * nbins: Size of above arrays
         *
         * Input:
         * dr2: transverse distance to spectra from particle (comoving kpc/h)
         * mass: mass of particle in absorping species (kg)
         * ppos: particle distance from box edge parallel to spectrum (comoving kpc/h)
         * pvel: particle velocity parallel to spectrum (physical km/s)
         * temp: particle temperature (K)
         * smooth: particle smoothing length (comoving kpc/h)
         */
        void add_particle(double * tau, double * colden, const int nbins, const double dr2, const double mass, const double ppos, const double pvel, const double temp, const double smooth);

        const double amumass, velfac, vbox;
};

#endif
