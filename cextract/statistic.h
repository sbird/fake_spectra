#ifndef SPB41_STATISTIC_H
#define SPB41_STATISTIC_H

#include <fftw3.h>
#include "global_vars.h"
/* These functions do the work*/
int powerspectrum(const int dims, double *field, double *power, fftw_plan* pl);
double mean_flux(double * tau, int nbins, double obs_flux, double tol);
void calc_power_spectra(double *flux_power, double *tau_H1,double scale, double tau_eff, int NumLos, int nbins);
void calc_pdf(double *flux_pdf, double *tau_H1,double scale, int NumLos, int nbins);
void smooth(double *in, double * out, int len,int slen);

void gaussian_smooth(double *in, double * out, int len,double sigma,double *kernel);
void smooth(double *in, double * out, int len,int slen);
double gaussian(double x, double sigma);
#endif
