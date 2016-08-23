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

/* Class extending LineAbsorption to efficiently loop
 * over all particles, using an index table to find those near each line
 */
class ParticleInterp: public LineAbsorption
{
    public:
        ParticleInterp(const int nbins_i, const double lambda, const double gamma, const double fosc, const double amumass, const double boxsize, const double velfac, const double atime, const double cofm[], const int axis[], const int NumLos, const int kernel):
        LineAbsorption(lambda, gamma, fosc, amumass, velfac, boxsize, atime, kernel),
        nbins(nbins_i),nlos(NumLos),
        sort_los_table(cofm, axis, NumLos, boxsize)
        {
        }

        /*Interpolate the particles to compute binned optical depth
         * at each bin along the line of sight. This is in redshift space,
         * including peculiar velocities.
         */
        void compute_tau(double * tau, const float Pos[], const float Vel[], const float Dens[], const float temp[], const float h[], const long long npart);
        /*Interpolate the particles to compute binned column density
         * at each bin along the line of sight. This is in physical space, not redshift space, ie,
         * peculiar velocities are ignored.
         */
        void compute_colden(double * colden, const float Pos[], const float Dens[], const float h[], const long long npart);

    private:
        const int nbins;
        const int nlos;
        IndexTable sort_los_table;
};

/*Convert the units of colden from atoms/(kpc/h)^2 to atoms (of species)/cm^2*/
void convert_colden_units(double * colden, const int nbins, const double h100, const double atime, const double amumass);

#endif
