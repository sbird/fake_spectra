#ifndef INDEX_TABLE_H
#define INDEX_TABLE_H

#include <map>
#include <vector>
#include <cmath>

#define RESO 0.1 // for arepo

/* The particles near each sightline, stored compressed: the entries for
 * line i are the range [offsets[i], offsets[i+1]) of part and dr2, in
 * increasing particle order.
 * A vector per line would need an allocation per line, and a map per line
 * an allocation per particle. This needs three, and can be filled in
 * parallel without a lock. */
struct NearParticles
{
    //One more than the number of lines: offsets[i+1] ends line i.
    std::vector<long long> offsets;
    //Particle index, and its squared distance from the line.
    std::vector<int> part;
    std::vector<double> dr2;

    //Number of particles near line i
    long long size(const int i) const
    {
        return offsets[i+1] - offsets[i];
    }

    int nlines() const
    {
        return offsets.size() - 1;
    }
};

class IndexTable
{
public:
  IndexTable(const double cofm[], const int axis[], const int NumLos_i, const double box);

  //Get a list of lines nearby a particle with coordinates xx, yy, zz and smoothing length hh.
  std::map<int,double> get_near_lines(const float pos[], const float hh);

  //Call fn(line index, squared distance) for each line near a particle.
  //Templated so that the callback is inlined into the search, which
  //means a caller need not build a container for each particle.
  template<class F> void for_each_near_line(const float pos[], const float hh, F fn)
  {
      if(index_table.size() > 0)
          each_nearby(pos[0], index_table, pos, hh, fn);
      if(index_table_xx.size() > 0)
          each_nearby(pos[1], index_table_xx, pos, hh, fn);
  }

  //Find the particles near each line.
  NearParticles get_near_particles(const float pos[], const float hh[], const long long npart);
  float * assign_cells(const int i, const NearParticles& nearby, const float pos[]);

  //Get the axis of a line
  inline int get_axis(const int iproc)
  {
      return axis[iproc];
  }

private:
  //Get the transverse distance from sightline iproc to position pos
  inline double calc_dr2(const double d1, const double d2)
  {
      /*    Distance to projection axis */
      double dr = fabs(d1);
      if(dr > 0.5*boxsize)
          dr = boxsize - dr; /* Keep dr between 0 and box/2 */
      double dr2 = dr*dr;
      dr = fabs(d2);
      if (dr > 0.5*boxsize)
          dr = boxsize - dr; /* between 0 and box/2 */
      dr2 += (dr*dr);
      return dr2;
  }

  inline bool second_close(const float second, const double lproj2, const float hh)
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

  //Visit the lines nearby a particle within an iterator range
  template<class F> void each_range(std::multimap<const double, const int>::const_iterator low, std::multimap<const double, const int>::const_iterator high, const float pos[], const float hh, const float first, F fn)
  {
      for(std::multimap<const double, const int>::const_iterator it = low; it != high; ++it)
      {
          const int iproc = it->second;
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
          if(second_close(second, lproj2, hh)){
              const double dr2 = calc_dr2(first - it->first, second - lproj2);
              if (dr2 <= hh*hh)
                  fn(iproc, dr2);
          }
      }
  }

  //Visit the lines nearby a particle from a particular index table
  template<class F> void each_nearby(float first, std::multimap<const double, const int>& sort_los, const float pos[], const float hh, F fn)
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
      //An iterator to the first element not less than ffm
      std::multimap<const double, const int>::const_iterator low = sort_los.lower_bound(ffm);
      //An iterator to the first element greater than ffp
      std::multimap<const double, const int>::const_iterator high = sort_los.lower_bound(ffp);
      //If periodic wrapping occurred, we want to go through zero
      if(ffm <= ffp) {
          each_range(low, high, pos, hh, first, fn);
      }
      else {
          each_range(sort_los.begin(), high, pos, hh, first, fn);
          each_range(low, sort_los.end(), pos, hh, first, fn);
      }
  }

  // The key is the position of the primary axis, which is xx for index_table and yy for index_table_xx.
  // The value is the index of this entry in cofm and axis.
  //index_table stores lines where axis = 2 or 3.
  std::multimap<const double, const int> index_table;
  //index_table_xx stores lines where axis = 1.
  std::multimap<const double, const int> index_table_xx;
  //Pointers to the original los table
  const double *cofm;
  const int *axis;
  const int NumLos;
  const double boxsize;
};

#endif
