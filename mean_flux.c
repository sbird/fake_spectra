/* Copyright (c) 2009, Simeon Bird <spb41@cam.ac.uk>
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

/*Calculates the scaling for the mean flux. 
 * tau is an array of optical depths, 
 * nbins is the total number of points in this array.
 * obs_flux is the value to scale the mean to.*/

#include <math.h>
#include <stdio.h>

double mean_flux(double * tau, double nbins, double obs_flux, double tol)
{
    double mean_flux=0;
    double tau_mean_flux=0;
    double scale=1;
    double newscale=1;
    double obs_flux_bins=obs_flux*nbins;
    double temp;
    int i;
    do{
         for(i=0; i< nbins; i++)
         {
             temp=exp(-scale*tau[i]);
             mean_flux+=temp;
             tau_mean_flux+=temp*tau[i];
         }
         newscale=scale+(mean_flux-obs_flux_bins)/tau_mean_flux;
         printf("iteration!\n");
    }while(fabs(newscale-scale) > tol*newscale);
    return newscale;
}
