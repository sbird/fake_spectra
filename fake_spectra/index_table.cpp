#include "index_table.h"
#include <stdio.h>
#include <cmath>
#include <cassert>
#ifdef _OPENMP
#include <omp.h>
#endif

//The thread that is running, and how many of them there could be.
//Compiling without OpenMP leaves one of each.
static inline int this_thread()
{
#ifdef _OPENMP
    return omp_get_thread_num();
#else
    return 0;
#endif
}

static inline int max_threads()
{
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}

static inline int cur_threads()
{
#ifdef _OPENMP
    return omp_get_num_threads();
#else
    return 1;
#endif
}

/* The two coordinates transverse to each (1-indexed) sightline axis. */
static const int PERP[4][2] = {{0,0},{1,2},{0,2},{0,1}};

// Bucket the lines of each axis onto their own mesh.
IndexTable::IndexTable(const double cofm_i[], const int axis_i[], const int NumLos_i, const double box):
  cofm(cofm_i), axis(axis_i), NumLos(NumLos_i), boxsize(box)
{
        std::vector<int> lines[3];
        for(int i=0;i<NumLos;i++){
            assert(axis[i] > 0 && axis[i] < 4);
            lines[axis[i]-1].push_back(i);
        }
        for(int ax = 0; ax < 3; ax++)
            mesh[ax].build(cofm, lines[ax], PERP[ax+1][0], PERP[ax+1][1], boxsize);
        return;
}

/*This function takes a particle list and returns, for each line, the list of
 * particles near it: near is defined as dx^2+dy^2 < h^2.
 * Each thread searches a contiguous chunk of the particles into its own
 * buffer, counting as it goes how many particles it found for each line.
 * Those counts give every thread a place to write in the output, so the
 * particles can be gathered without a lock, and in one pass over them.
 * Each thread searches its chunk in order, and the chunks are handed out
 * in order, so the particles of a line come out sorted by particle index,
 * as they were when each line kept a std::map. */
NearParticles IndexTable::get_near_particles(const float pos[], const float hh[], const long long npart)
{
    const int nthread = max_threads();
    //Particle indices are 64 bit, line indices are bounded by NumLos.
    std::vector<std::vector<long long> > tpart(nthread);
    std::vector<std::vector<int> > tline(nthread);
    std::vector<std::vector<double> > tdr2(nthread);
    //How many particles each thread found for each line.
    std::vector<long long> counts((long long)nthread*NumLos, 0);
    #pragma omp parallel
    {
        const int tid = this_thread();
        std::vector<long long>& mypart = tpart[tid];
        std::vector<int>& myline = tline[tid];
        std::vector<double>& mydr2 = tdr2[tid];
        long long * mycount = &counts[(long long)tid*NumLos];
        /*Most particles reaching here are near one line, so this is about the
         * right size, and saves regrowing the buffers as they fill.*/
        const size_t guess = npart/nthread + 16;
        mypart.reserve(guess);
        myline.reserve(guess);
        mydr2.reserve(guess);
        /* Static schedule hands the chunks out round robin in order of the thread number.
         * So thread t searches the t'th chunk, which is what puts the particles of a line in order.*/
        const int nth = cur_threads();
        const long long chunk = (npart + nth - 1)/nth;
        #pragma omp for schedule(static, chunk)
        for(long long i = 0; i < npart; i++){
            for_each_near_line(&(pos[3*i]), hh[i], [&](const int iproc, const double dr2){
                    mypart.push_back(i);
                    myline.push_back(iproc);
                    mydr2.push_back(dr2);
                    mycount[iproc]++;
                });
        }
    }
    //Where each thread starts writing for each line: line major and thread
    //minor, so that the chunks of a line end up in the order searched.
    NearParticles nearby;
    nearby.offsets.resize(NumLos+1);
    std::vector<long long> start((long long)nthread*NumLos);
    long long total = 0;
    for(int il = 0; il < NumLos; il++){
        nearby.offsets[il] = total;
        for(int it = 0; it < nthread; it++){
            start[(long long)it*NumLos + il] = total;
            total += counts[(long long)it*NumLos + il];
        }
    }
    nearby.offsets[NumLos] = total;
    nearby.part.resize(total);
    nearby.dr2.resize(total);
    //One iteration per buffer, however many threads there turn out to be.
    #pragma omp parallel for schedule(static, 1)
    for(int it = 0; it < nthread; it++){
        long long * where = &start[(long long)it*NumLos];
        const size_t nent = tpart[it].size();
        for(size_t k = 0; k < nent; k++){
            const long long ii = where[tline[it][k]]++;
            nearby.part[ii] = tpart[it][k];
            nearby.dr2[ii] = tdr2[it][k];
        }
    }
    return nearby;
}

float * IndexTable::assign_cells(const int line_i, const NearParticles& nearby, const float pos[])
{
    const long long first = nearby.offsets[line_i];
    const long long Ncells = nearby.size(line_i);
    // printf("assigning parts of line %d to %d cells...\n", line_i, Ncells);
    float * arr2 = new float [2*Ncells];
    //Nothing is near this line, so there is nothing to assign the
    //grid points to. The caller's particle loop is empty as well.
    if(Ncells == 0)
        return arr2;
    // initialize
    for(long long i = 0; i < 2*Ncells; ++i)
        arr2[i] = 3*boxsize;

    // divide each sightline into an array. grid size = RESO ckpc/h
    const int N = int(boxsize/RESO);
    const double reso = boxsize/N;

    const int axis_i = axis[line_i];
    const double yp = cofm[3*line_i+axis_i%3], zp = cofm[3*line_i+(axis_i+1)%3];
    // loop over all grid points along the sightline
    for(int i = 0; i < N; ++i)
    {
        // the default choice is the sightline points in the x direction
        double xp = (i+0.5)*reso;
        // find the particle index that this point along the sightline belongs to
        double min_dist = boxsize;
        long long min_ind = 0;
        for(long long ind = 0; ind < Ncells; ++ind){
            const long long ipart = nearby.part[first+ind];
            double dx, dy, dz;
            // take into account periodicity
            dx = fabs(pos[3*ipart+axis_i-1]-xp);
            if(dx > boxsize/2.) dx = boxsize - dx;
            dy = fabs(pos[3*ipart+axis_i%3]-yp);
            if(dy > boxsize/2.) dy = boxsize - dy;
            dz = fabs(pos[3*ipart+(axis_i+1)%3]-zp);
            if(dz > boxsize/2.) dz = boxsize - dz;
            double dist = sqrt(dx*dx + dy*dy + dz*dz);
            if(dist < min_dist){
                min_dist = dist;
                min_ind = ind;
            }
        }
        //This changes sign in the special
        //case where we have wrapped around the box.
        //We can break here because we know all the remaining
        //grid points will be close to this particle.
        if(arr2[2*min_ind] < reso && xp > boxsize/2.+0.5*reso)
        {
            arr2[2*min_ind] = xp;
            arr2[2*min_ind+1] += boxsize;
            break;
        }
        //This asserts that we are advancing the maximal
        //cell only one grid point.
        if(xp > 1.5*reso + arr2[2*min_ind+1])
        {
            printf("Advanced pointer more than expected for cell: %d left=%f, right=%f\n", i, arr2[2*min_ind], arr2[2*min_ind+1]);
            exit(1);
        }
        // Assign lowermost grid point if this is the first place
        // this particle is close to.
        if(arr2[2*min_ind] > 2*boxsize)
            arr2[2*min_ind] = xp;
        //Assign uppermost grid point
        arr2[2*min_ind+1] = xp;
    }

    for(long long i = 0; i < Ncells; ++i){
        arr2[2*i] -= 0.5*reso;
        arr2[2*i+1] += 0.5*reso;
    }

    return arr2;
}
