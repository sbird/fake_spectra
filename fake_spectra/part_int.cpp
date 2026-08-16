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

#include "part_int.h"
#include "singleabs.h"
//For NULL
#include <cstddef>

void ParticleInterp::compute_tau(double tau[], const double Pos[], const double Vel[], const double Dens[], const double temp[], const double h[], const long long npart)
{
    const NearParticles nearby_array = sort_los_table.get_near_particles(Pos, h, npart);
    const long long nlines = nearby_array.nlines();
    #pragma omp parallel for
    for(long long i = 0; i < nlines; ++i)
    {
        double * arr2 = NULL;
        if(kernel == VORONOI_MESH) arr2 = sort_los_table.assign_cells(i, nearby_array, Pos);
        const int axis = sort_los_table.get_axis(i);
        double * tau_loc = &tau[i*nbins];
        //List of particles near this los
        //Loop over them
        const long long first = nearby_array.offsets[i];
        const long long nnear = nearby_array.size(i);
        for(long long ind = 0; ind < nnear; ++ind)
        {
          const long long ipart = nearby_array.part[first+ind];
          const double dr2 = nearby_array.dr2[first+ind];
          //Particle position parallel to axis
          const double ppos = Pos[3*ipart+axis-1];
          const double pvel = Vel[3*ipart+axis-1];
          if(kernel == VORONOI_MESH)
              add_tau_particle(tau_loc, nbins, arr2[2*ind], Dens[ipart], ppos, pvel, temp[ipart], arr2[2*ind+1]);
          else
              add_tau_particle(tau_loc, nbins, dr2, Dens[ipart], ppos, pvel, temp[ipart], h[ipart]);
        }  /*Loop over list of particles near LOS*/
        if(kernel == VORONOI_MESH) delete [] arr2;
    } /* Loop over LOS*/
    return;
}

void ParticleInterp::compute_colden(double colden[], const double Pos[], const double Dens[], const double h[], const long long npart)
{
    const NearParticles nearby_array = sort_los_table.get_near_particles(Pos, h, npart);
    const long long nlines = nearby_array.nlines();
    #pragma omp parallel for
    for(long long i = 0; i < nlines; ++i)
    {
        double * arr2 = NULL;
        if(kernel == VORONOI_MESH) arr2 = sort_los_table.assign_cells(i, nearby_array, Pos);
        const int axis = sort_los_table.get_axis(i);
        double * colden_loc = &colden[i*nbins];
        //List of particles near this los
        //Loop over them
        const long long first = nearby_array.offsets[i];
        const long long nnear = nearby_array.size(i);
        for(long long ind = 0; ind < nnear; ++ind)
        {
          const long long ipart = nearby_array.part[first+ind];
          const double dr2 = nearby_array.dr2[first+ind];
          //Particle position parallel to axis
          const double ppos = Pos[3*ipart+axis-1];
          //Don't need temp if no tau
          if(kernel == VORONOI_MESH)
              add_colden_particle(colden_loc, nbins, arr2[2*ind], Dens[ipart], ppos, arr2[2*ind+1]);
          else
              add_colden_particle(colden_loc, nbins, dr2, Dens[ipart], ppos, h[ipart]);
        }  /*Loop over list of particles near LOS*/
        if(kernel == VORONOI_MESH) delete [] arr2;
    } /* Loop over LOS*/
    return;
}
