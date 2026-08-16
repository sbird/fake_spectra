#ifndef INDEX_TABLE_H
#define INDEX_TABLE_H

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
    //The index is 64 bit: a snapshot may have more than 2^31 particles.
    std::vector<long long> part;
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

/* The sightlines with one particular axis, bucketed onto a uniform grid in
 * the two coordinates transverse to that axis.
 * The lines of a cell are contiguous, and each one's transverse coordinates
 * are stored beside its index, so a search reads only the cells it needs,
 * in order, and never reads back into the cofm table.
 * The grid holds about one line per cell, which means the number of lines a
 * search looks at is set by how densely the lines are spread and by the
 * smoothing length, and not by how many lines there are in total.
 * Indexing both transverse coordinates is the point: a table sorted on one
 * of them has to look at every line within a smoothing length in that
 * coordinate, which is a number that grows with the number of lines. */
class LineMesh
{
public:
    LineMesh(): ncell(0), boxsize(0), invcell(0), p1(0), p2(0) {}

    /*The lines given are those with one axis; p1 and p2 are the two
     * coordinates transverse to it.*/
    void build(const double cofm[], const std::vector<int>& lines, const int p1_i, const int p2_i, const double box)
    {
        p1 = p1_i;
        p2 = p2_i;
        boxsize = box;
        const size_t nlines = lines.size();
        ncell = (int) sqrt((double)nlines);
        if(ncell < 1)
            ncell = 1;
        invcell = ncell/boxsize;
        cellstart.assign((size_t)ncell*ncell+1, 0);
        c1.resize(nlines);
        c2.resize(nlines);
        lineid.resize(nlines);
        //Counting sort of the lines into their cells.
        std::vector<int> cell(nlines);
        for(size_t i = 0; i < nlines; i++) {
            cell[i] = which_cell(cofm[3*lines[i]+p1], cofm[3*lines[i]+p2]);
            cellstart[cell[i]+1]++;
        }
        for(size_t i = 1; i < cellstart.size(); i++)
            cellstart[i] += cellstart[i-1];
        std::vector<int> where(cellstart.begin(), cellstart.end()-1);
        for(size_t i = 0; i < nlines; i++) {
            const int k = where[cell[i]]++;
            //Kept in double, as they are in cofm: rounding the line position
            //to float would move dr2 in the last few digits.
            c1[k] = cofm[3*lines[i]+p1];
            c2[k] = cofm[3*lines[i]+p2];
            lineid[k] = lines[i];
        }
    }

    bool empty() const
    {
        return lineid.empty();
    }

    //Call fn(line index, squared distance) for each line within hh of the particle.
    template<class F> void query(const double pos[], const double hh, F fn) const
    {
        const double x = pos[p1], y = pos[p2];
        int i0 = ifloor((x-hh)*invcell), i1 = ifloor((x+hh)*invcell);
        int j0 = ifloor((y-hh)*invcell), j1 = ifloor((y+hh)*invcell);
        //A particle reaching more than a box covers every cell, and wrapping
        //the cell index would visit some of them more than once.
        if(i1 - i0 + 1 >= ncell) { i0 = 0; i1 = ncell-1; }
        if(j1 - j0 + 1 >= ncell) { j0 = 0; j1 = ncell-1; }
        const double hh2 = hh*hh;
        for(int i = i0; i <= i1; i++) {
            const int ii = wrap(i);
            for(int j = j0; j <= j1; j++) {
                const int jj = wrap(j);
                const int c = ii*ncell + jj;
                const int end = cellstart[c+1];
                for(int k = cellstart[c]; k < end; k++) {
                    const double dr2 = calc_dr2(x - c1[k], y - c2[k]);
                    if(dr2 <= hh2)
                        fn(lineid[k], dr2);
                }
            }
        }
    }

private:
    inline int ifloor(const double x) const
    {
        const int i = (int) x;
        return (x < i) ? i-1 : i;
    }

    //The cell range searched extends at most one box either way, so the
    //index is at most one box outside the grid.
    inline int wrap(const int i) const
    {
        if(i < 0)
            return i + ncell;
        if(i >= ncell)
            return i - ncell;
        return i;
    }

    inline int which_cell(const double x, const double y) const
    {
        int i = (int)(x*invcell);
        int j = (int)(y*invcell);
        if(i < 0) i = 0;
        if(i >= ncell) i = ncell-1;
        if(j < 0) j = 0;
        if(j >= ncell) j = ncell-1;
        return i*ncell + j;
    }

    //Get the transverse distance from a sightline to a position
    inline double calc_dr2(const double d1, const double d2) const
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

    int ncell;
    double boxsize, invcell;
    //Which two coordinates of a position this mesh is indexed on.
    int p1, p2;
    //Where each cell starts in the three arrays below: ncell*ncell+1 entries.
    std::vector<int> cellstart;
    //The transverse coordinates of each line and its index in cofm and axis,
    //in cell order.
    std::vector<double> c1, c2;
    std::vector<int> lineid;
};

class IndexTable
{
public:
  IndexTable(const double cofm[], const int axis[], const int NumLos_i, const double box);

  //Find the particles near each line.
  NearParticles get_near_particles(const double pos[], const double hh[], const long long npart);
  double * assign_cells(const int i, const NearParticles& nearby, const double pos[]);

  //Get the axis of a line
  inline int get_axis(const int iproc)
  {
      return axis[iproc];
  }

private:
  //Call fn(line index, squared distance) for each line near a particle.
  //Templated so that the callback is inlined into the search, which
  //means the caller need not build a container for each particle.
  template<class F> void for_each_near_line(const double pos[], const double hh, F fn)
  {
      for(int ax = 0; ax < 3; ax++)
          if(!mesh[ax].empty())
              mesh[ax].query(pos, hh, fn);
  }

  //One mesh per sightline axis, so that every line in a mesh is indexed on
  //the same two coordinates and the search needs no per-line branch.
  LineMesh mesh[3];
  //Pointers to the original los table
  const double *cofm;
  const int *axis;
  const int NumLos;
  const double boxsize;
};

#endif
