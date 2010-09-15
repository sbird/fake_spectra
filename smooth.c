#include "statistic.h"
#include <math.h>

double gaussian(double x, double sigma){
        return exp(-x*x/(2*sigma*sigma))/(sqrt(2*M_PI)*sigma);
}

int pstv(int a){
        return a>0 ? a: -a;
}
void gaussian_smooth(double *in, double * out, int len,double sigma,double * kernel)
{
    int i,j;
    int off=3*sigma;
    /*Smooth the middle*/
    for(i=0;i<len; i++){
       out[i]=0;
       for(j=-off; j<=off; j++){
                int ind=(i+j)%len;
                if(ind<0) ind+=len;
                out[i]+=in[ind]*kernel[pstv(j)];
       }
    }
    return;
}


void smooth(double *in, double * out, int len,int slen)
{
    int i,j;
    int off=(slen-1)/2;
    /*For the start*/
    for(i=0;i<off;i++){
        out[i]=0;
        /*Pad out by repeating first value*/
        for(j=-off;j<-i;j++)
            out[i]+=in[0];
        for(j=-i;j<=off;j++)
            out[i]+=in[i+j];
    }
    /*Smooth the middle*/
    for(i=off;i<len-off; i++){
       out[i]=0;
       for(j=-off; j<=off; j++)
                out[i]+=in[i+j];
    }
    /*For the end*/
    for(i=len-off;i<len;i++){
        out[i]=0;
        for(j=-off;j<len-i;j++)
            out[i]+=in[i+j];
        /*Pad out by repeating last value*/
        for(j=len-i;j<=-i;j++)
            out[i]+=in[len-1];
    }
    for(i=0;i<len;i++)
            out[i]/=slen;
    return;
}

