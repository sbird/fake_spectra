#include "index_table.h"
#include <cmath>
#include <cassert>

// Construct the index tables as maps, which are automatically sorted.
IndexTable::IndexTable(const double cofm_i[], const int axis_i[], const int NumLos_i, const double box):
  cofm(cofm_i), axis(axis_i), NumLos(NumLos_i), boxsize(box)
{
        for(int i=0;i<NumLos;i++){
            if(axis[i] == 1)
                index_table_xx.insert(std::pair<const double, const int>(cofm[3*i+1], i));
            else
                index_table.insert(std::pair<const double, const int>(cofm[3*i], i));
        }
        assert(index_table_xx.size() + index_table.size() == (unsigned int) NumLos);
        return;
}

/*Returns a std::map of lines close to the coordinates xx, yy, zz.
 * the key is the line index, and the value is the distance from the two axes not projected along*/
void IndexTable::get_nearby_from_range(std::multimap<const double, const int>::const_iterator low, std::multimap<const double, const int>::const_iterator high, std::map<int, double>& nearby, const float pos[], const float hh, const float first)
{
      for(std::multimap<const double, const int>::const_iterator it = low; it != high;++it)
      {
          const int iproc=it->second;
          /*If close in the second coord, save line*/
          /*Load a sightline from the table.*/
          const int iaxis = axis[iproc];

          float second;
          double lproj2;
          if (iaxis == 3){
            second = pos[1];
            lproj2 = cofm[3*iproc+1];
          }
          else{
            second = pos[2];
            lproj2 = cofm[3*iproc+2];
          }
          const double lproj = it->first;

          if(second_close(second, lproj2, hh)){
              double dr2 = calc_dr2(first-lproj, second-lproj2);
              if (dr2 <= hh*hh){
                      nearby[iproc]=dr2;
              }
          }
      }
}

inline bool IndexTable::second_close(const float second, const double lproj2, const float hh)
{
    /* Now check that xx-hh < proj < xx +  hh */
    float ffp=second+hh;
    //Periodic wrap
    if(ffp > boxsize)
        if(lproj2 < ffp - boxsize)
            return true;
    float ffm=second-hh;
    if(ffm < 0)
        if(lproj2 > ffm + boxsize)
            return true;
    if (lproj2 > ffm && lproj2 < ffp)
        return true;
    else
        return false;
}

inline double IndexTable::calc_dr2(const double d1, const double d2)
{
    double dr, dr2;
    /*    Distance to projection axis */
    dr = fabs(d1);

    if(dr > 0.5*boxsize)
            dr = boxsize - dr; /* Keep dr between 0 and box/2 */

    dr2 = dr*dr;

    dr = fabs(d2);
    if (dr > 0.5*boxsize)
      dr = boxsize - dr; /* between 0 and box/2 */

    dr2 += (dr*dr);
    return dr2;
}

void IndexTable::get_nearby(float first, std::multimap<const double, const int>& sort_los, std::map<int, double>& nearby, const float pos[], const float hh)
{
      /*Now find the elements where dr < 2 hh, wrapping with respect to boxsize*/
      /* First find highest index where xx + 2 hh > priax */
      float ffp=first+hh;
      if(ffp > boxsize)
         ffp-=boxsize;
      /* Now find lowest index in what remains where xx - 2 hh < priax */
      float ffm=first-hh;
      if(ffm < 0)
         ffm+=boxsize;
      std::multimap<const double, const int>::const_iterator low,high;
      //An iterator to the first element not less than ffm
      low=sort_los.lower_bound(ffm);
      //An iterator to the first element greater than ffp
      high=sort_los.lower_bound(ffp);
      //If periodic wrapping occurred, we want to go through zero
      if(ffm <= ffp) {
        get_nearby_from_range(low, high, nearby, pos, hh, first);
      }
      else {
        get_nearby_from_range(sort_los.begin(), high, nearby, pos, hh, first);
        get_nearby_from_range(low, sort_los.end(), nearby, pos, hh, first);
      }
}

/*This function takes a particle position and returns a list of the indices of lines near it in index_nr_lines
 * Near is defined as: dx^2+dy^2 < 4h^2 */
std::map<int,double> IndexTable::get_near_lines(const float pos[],const float hh)
{
      std::map<int, double> nearby;
      if(index_table.size() > 0){
        get_nearby(pos[0],index_table, nearby, pos, hh);
      }
      if(index_table_xx.size() > 0){
        get_nearby(pos[1],index_table_xx, nearby, pos, hh);
      }
      return nearby;
}

//Find a list of particles near each line
std::valarray< std::map<int, double> > IndexTable::get_near_particles(const float pos[], const float hh[], const long long npart)
{
    //List of lines. Each element contains a list of particles and their distances to the line.
    std::valarray< std::map<int, double> > line_near_parts(NumLos);
    //find lists
    #pragma omp parallel for
    for(long long i=0; i < npart; i++){
        //Get list of lines near this particle
	    std::map<int, double> nearby=get_near_lines(&(pos[3*i]),hh[i]);

        if(nearby.size()){
            #pragma omp critical
            {
                //Insert the particles into the list of particle lists
                for(std::map <int, double>::const_iterator it = nearby.begin(); it != nearby.end(); ++it)
                       line_near_parts[it->first].insert(std::pair<int, double>(i, it->second));
            }
        }
    }
    return line_near_parts;
}
