#include "statistic.h"
#include "parameters.h"

void calc_power_spectra(double *flux_power, double *tau_H1,double scale,double tau_eff,int NumLos)
{
    int j,i;
    double flux_power_local[(NBINS+1)/2];
    double flux_power_sm[(NBINS+1)/2];
    double flux_H1_local[NBINS];
    fftw_plan pl=rfftw_create_plan(NBINS,FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE);
    for(j=0; j<(NBINS+1)/2;j++)
      flux_power[j]=0;
    /*Perform the scaling*/
    for(j=0; j<NBINS*NumLos; j++)
      tau_H1[j]*=scale;
    for(j=1;j<NumLos-1;j++)
    {
        /* Calculate flux and flux power spectrum */
        for(i=0; i<NBINS; i++)
           flux_H1_local[i]=exp(-tau_H1[j*NBINS+i])/exp(-tau_eff)-1.;
        powerspectrum(NBINS, flux_H1_local, flux_power_local, pl);
        /*Smooth*/
        flux_power_sm[0]=flux_power_local[0];
        smooth(flux_power_local+1,flux_power_sm+1,(NBINS+1)/2-1,3);
        /*Write powerspectrum*/
        for(i=0; i<(NBINS+1)/2;i++)
            flux_power[i]+=flux_power_sm[i];
    }
    /*Average the power spectrum*/
    for(j=0; j<(NBINS+1)/2;j++)
        flux_power[j]/=(NumLos-2);
    fftw_destroy_plan(pl);
    return;
}


