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

void ParticleInterp::compute_tau(double tau[], const float Pos[], const float Vel[], const float Dens[], const float temp[], const float h[], const long long npart)
{
    const std::valarray< std::map<int, double> > nearby_array = sort_los_table.get_near_particles(Pos, h, npart);
    //Use a plain int as not sure openmp can handle iterators efficiently.
    const unsigned int nlines = nearby_array.size();
    #pragma omp parallel for
    for(unsigned int i = 0; i < nlines; ++i)
    {
        float * arr2 = NULL;
        if(kernel == VORONOI_MESH) arr2 = sort_los_table.assign_cells(i, nearby_array, Pos);
        const int axis = sort_los_table.get_axis(i);
        double * tau_loc = &tau[i*nbins];
        //List of particles near this los
        //Loop over them
        int ind = 0;
        for(std::map<int, double>::const_iterator it = nearby_array[i].begin(); it != nearby_array[i].end(); ++it)
        {
          const int ipart = it->first;
          const double dr2 = it->second;
          //Particle position parallel to axis
          const float ppos = Pos[3*ipart+axis-1];
          const float pvel = Vel[3*ipart+axis-1];
          if(kernel == VORONOI_MESH)
              add_tau_particle(tau_loc, nbins, arr2[2*ind], Dens[ipart], ppos, pvel, temp[ipart], arr2[2*ind+1]);
          else
              add_tau_particle(tau_loc, nbins, dr2, Dens[ipart], ppos, pvel, temp[ipart], h[ipart]);
          ind++;
        }  /*Loop over list of particles near LOS*/
        if(kernel == VORONOI_MESH) delete [] arr2;
    } /* Loop over LOS*/
    return;
}

void ParticleInterp::compute_colden(double colden[], const float Pos[], const float Dens[], const float h[], const long long npart)
{
    const std::valarray< std::map<int, double> > nearby_array = sort_los_table.get_near_particles(Pos, h, npart);
    //Use a plain int as not sure openmp can handle iterators efficiently.
    const unsigned int nlines = nearby_array.size();
    #pragma omp parallel for
    for(unsigned int i = 0; i < nlines; ++i)
    {
        float * arr2 = NULL;
        if(kernel == VORONOI_MESH) arr2 = sort_los_table.assign_cells(i, nearby_array, Pos);
        const int axis = sort_los_table.get_axis(i);
        double * colden_loc = &colden[i*nbins];
        //List of particles near this los
        //Loop over them
        int ind = 0;
        for(std::map<int, double>::const_iterator it = nearby_array[i].begin(); it != nearby_array[i].end(); ++it)
        {
          const int ipart = it->first;
          const double dr2 = it->second;
          //Particle position parallel to axis
          const float ppos = Pos[3*ipart+axis-1];
          //Don't need temp if no tau
          if(kernel == VORONOI_MESH)
              add_colden_particle(colden_loc, nbins, arr2[2*ind], Dens[ipart], ppos, arr2[2*ind+1]);
          else
              add_colden_particle(colden_loc, nbins, dr2, Dens[ipart], ppos, h[ipart]);
          ind++;
        }  /*Loop over list of particles near LOS*/
        if(kernel == VORONOI_MESH) delete [] arr2;
    } /* Loop over LOS*/
    return;
}
