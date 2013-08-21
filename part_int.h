#ifndef PART_INT_H
#define PART_INT_H

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

#include "index_table.h"
#include "absorption.h"

/* Class extending ComputeLineAbsorption to efficiently loop
 * over all particles, using an index table to find those near each line
 */
class ParticleInterp: public ComputeLineAbsorption
{
    public:
        ParticleInterp(double * tau_i, double * colden_i, const int nbins_i, const double lambda, const double gamma, const double fosc, const double amumass, const double boxsize, const double velfac, const los *los_table_i, const int NumLos):
        ComputeLineAbsorption(lambda, gamma, fosc, amumass, velfac, boxsize),
        tau(tau_i), colden(colden_i),nbins(nbins_i),nlos(NumLos),
        sort_los_table(los_table_i, NumLos, boxsize)
        {
        }

        /*Interpolate the particles in the given arrays onto
         * the spectra pointed at by tau and colden
         * */
        void do_work(const float Pos[], const float Vel[], const float Mass[], const float temp[], const float h[], const long long npart);

    private:
        double *tau, *colden;
        const int nbins;
        const int nlos;
        IndexTable sort_los_table;
};

/*Convert the units of colden from atoms/(kpc/h)^2 to atoms/cm^2*/
void convert_colden_units(double * colden, const int nbins, const double h100, const double atime);

#endif
