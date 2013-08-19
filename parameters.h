#ifndef PARAMETERS_H
#define PARAMETERS_H

/* Snapshot information */
#define PARTTYPE 0/* Particle type required */
#define FILENUMBER 1 /* Number of sub-files */

/* Spectrum data set size     */ 
#define NBINS 1024 /* number of pixels */

/* Model parameters outwith header */
#define XH 0.76  /* hydrogen fraction by mass */
/*The value from 0711.1862 is (0.0023±0.0007) (1+z)^(3.65±0.21)*/
#define TAU_EFF 0.0023*pow(1.0+redshift,3.65)

/* Some useful numbers */
#define  GAMMA (5.0/3.0)

/* Physical constants, SI units */
#define  GRAVITY     6.67428e-11
#define  KPC 3.08568025e19
#define  BOLTZMANN   1.3806504e-23  /* m2 kg s-2 K-1 */
#define  PROTONMASS  1.66053886e-27 /* 1 a.m.u */
#define  SOLAR_MASS 1.98892e30
#define  LIGHT           2.99792458e8 /*in km/s*/
#define  SIGMA_T 6.652458558e-29 /* Thompson cross-section in m^2*/

/* Atomic data (from VPFIT) */
#define  LAMBDA_LYA_H1 1215.6701e-10
#define  LAMBDA_LYA_HE2 303.7822e-10
#define  FOSC_LYA 0.416400
#define  HMASS 1.00794   /* Hydrogen mass in a.m.u. */
#define  HEMASS 4.002602 /* Helium-4 mass in a.m.u. */
#define  GAMMA_LYA_H1 6.265e8
#define  GAMMA_LYA_HE2 6.27e8

// convert energy/unit mass to J kg^-1
#define  ESCALE (1.0e6)

#endif
