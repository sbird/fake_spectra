#include <math.h>
#include "statistic.h"
#include <math.h>

/* #define SIGMA 2.5 */

void calc_power_spectra(double *flux_power, double *tau_H1,double scale,double tau_eff,int NumLos)
{
    int j,i;
    double flux_power_local[(NBINS+1)/2];
/*     double flux_H1_sm[NBINS]; */
    double flux_power_sm[(NBINS+1)/2];
    double flux_H1_local[NBINS+2];
/*     double kernel[(int)ceil(3*SIGMA+1)]; */
    /*Calculate kernel table*/
/*     for(i=0;i<=3*SIGMA;i++){ */
/*            kernel[i]=gaussian(i,SIGMA); */
/*     } */
    fftw_plan pl=fftw_plan_dft_r2c_1d(NBINS,flux_H1_local,(fftw_complex *)flux_H1_local, FFTW_ESTIMATE);
    for(j=0; j<(NBINS+1)/2;j++)
      flux_power[j]=0;
    /*Perform the scaling*/
    for(j=0; j<NBINS*NumLos; j++)
      tau_H1[j]*=scale;
    for(j=0;j<NumLos;j++)
    {
        /* Calculate flux and flux power spectrum */
        for(i=0; i<NBINS; i++)
           flux_H1_local[i]=exp(-tau_H1[j*NBINS+i])/exp(-tau_eff)-1.;
        /*Smooth*/
/*         gaussian_smooth(flux_H1_local,flux_H1_sm,NBINS,SIGMA,kernel); */
        powerspectrum(NBINS, flux_H1_local, flux_power_local, &pl);
        flux_power_sm[0]=flux_power_local[0];
        smooth(flux_power_local+1,flux_power_sm+1,(NBINS+1)/2-1,3);
        /*Write powerspectrum*/
        for(i=0; i<(NBINS+1)/2;i++)
            flux_power[i]+=flux_power_sm[i];
    }
    /*Average the power spectrum*/
    for(j=0; j<(NBINS+1)/2;j++)
        flux_power[j]/=(NumLos);//*exp(-pow(j*6.25*SIGMA/(NBINS*1.0),2)));
    fftw_destroy_plan(pl);
    return;
}


