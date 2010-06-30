#include <drfftw.h>
/* These functions do the work*/
int powerspectrum(const int dims, double *field, double *power, fftw_plan pl);
double mean_flux(double * tau, int nbins, double obs_flux, double tol);
void calc_power_spectra(double *flux_power, double *tau_H1,double scale, double tau_eff, int NumLos);
void calc_pdf(double *flux_pdf, double *tau_H1,double scale, int NumLos);
