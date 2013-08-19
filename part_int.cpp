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


void ParticleInterp::do_work(const float Pos[], const float Vel[], const float Mass[], const float temp[], const float h[], const int npart)
{
    for(int i=0;i<npart;i++)
    {
      std::map<int, double> nearby=sort_los_table.get_near_lines(&(Pos[3*i]),h[i]*0.5);
      for(std::map<int, double>::iterator it = nearby.begin(); it != nearby.end(); ++it)
      {
          const int iproc=it->first;
          const double dr2=it->second;
          /*Load a sightline from the table.*/
          const int iaxis = los_table[iproc].axis;
          //Particle position parallel to axis
          const float ppos = Pos[3*i+iaxis-1];
          const float pvel = Vel[3*i+iaxis-1];
          double * tau_loc = (tau ? &tau[iproc*nbins] : NULL);
          add_particle(tau_loc, &colden[iproc*nbins], nbins, dr2, Mass[i], ppos, pvel, temp[i], h[i]);
	  }  /*Loop over LOS*/
    } /* Loop over particles*/
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
