#include "headers.h"


/* Snapshot information */
#define PARTTYPE 0/* Particle type required */
#define FILENUMBER 1 /* Number of sub-files */


/* Optional flags, set to either 1 (on) or 0 (off) */
#define PERIODIC 1 /* Periodic spectra */
#define VOIGT 1    /* Voigt profiles vs. Gaussian profiles */
#define PECVEL 1   /* Peculiar velocities */

/* Spectrum data set size     */ 
#define NBINS 1024 /* number of pixels */
#define NUMLOS 20 /* number of spectra */

/* Model parameters outwith header */
#define XH 0.76  /* hydrogen fraction by mass */
#define OMEGAB 0.0463 /* baryon fraction */

/* Some useful numbers */
#define  PI    3.14159265358979323846
#define  GAMMA (5.0/3.0)

/* Physical constants, SI units */
#define  GRAVITY     6.67428e-11
#define  BOLTZMANN   1.3806504e-23
#define  C           2.99792458e8
#define  PROTONMASS  1.66053886e-27 /* 1 a.m.u */
#define  MPC 3.08568025e22
#define  KPC 3.08568025e19
#define  SIGMA_T 6.652458558e-29
#define  SOLAR_MASS 1.98892e30

/* Atomic data (from VPFIT) */
#define  LAMBDA_LYA_H1 1215.6701e-10
#define  LAMBDA_LYA_HE2 303.7822e-10
#define  FOSC_LYA 0.416400
#define  HMASS 1.00794   /* Hydrogen mass in a.m.u. */
#define  HEMASS 4.002602 /* Helium-4 mass in a.m.u. */
#define  GAMMA_LYA_H1 6.265e8
#define  GAMMA_LYA_HE2 6.27e8
