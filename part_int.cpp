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
//For NULL
#include <cstddef>

void ParticleInterp::do_work(const float Pos[], const float Vel[], const float Mass[], const float temp[], const float h[], const long long npart)
{
    const std::valarray< std::vector<std::pair<int, double> > > nearby_array = sort_los_table.get_near_particles(Pos, h, npart);
    //Use a plain int as not sure openmp can handle iterators efficiently.
    //Parallelizing this bit makes little difference: most of the time is spent in get_near_particles
    const unsigned int nlines = nearby_array.size();
    #pragma omp parallel for
    for(unsigned int i = 0; i < nlines; ++i)
    {
        const int axis = sort_los_table.get_axis(i);
        double * tau_loc = (tau ? &tau[i*nbins] : NULL);
        double * colden_loc = &colden[i*nbins];
        //List of particles near this los
        //Loop over them
        for(std::vector<std::pair<int, double> >::const_iterator it = nearby_array[i].begin(); it != nearby_array[i].end(); ++it)
        {
          const int ipart = it->first;
          const double dr2 = it->second;
          //Particle position parallel to axis
          const float ppos = Pos[3*i+axis-1];
          const float pvel = Vel[3*i+axis-1];
          add_particle(tau_loc, colden_loc, nbins, dr2, Mass[ipart], ppos, pvel, temp[ipart], h[ipart]);
	    }  /*Loop over list of particles near LOS*/
    } /* Loop over LOS*/
    return;
}

/*Convert the units of colden from atoms/(kpc/h)^2 to atoms/cm^2*/
void convert_colden_units(double * colden, const int nbins, const double h100, const double atime)
{
  /* Conversion factors from internal units */
  const double rscale = (KPC*atime)/h100;
  for(int i = 0;i<nbins;i++)
     colden[i]*=pow(rscale,-2);
  return;
}
