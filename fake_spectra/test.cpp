/* Copyright (c) 2013, Simeon Bird <spb@ias.edu>
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
#define BOOST_TEST_DYN_LINK

/** \file
 * Test suite using boost::test*/
#define BOOST_TEST_MODULE FluxExtract
#include "absorption.h"
#include "singleabs.h"
#include "index_table.h"
#include <cstdio>
#include <numeric>
#include <math.h>
#include <set>
#include <boost/test/unit_test.hpp>
#include <boost/test/test_tools.hpp>

#define FLOATS_NEAR_TO(x,y) \
        BOOST_CHECK_MESSAGE( !isinf((x)) && !isinf((y))  && !isnan((x)) && !isnan((y)) \
                && fabs((x) - (y)) <= std::max<double>(fabs(x),fabs(y))/1e5 ,(x)<<" is not close to "<<(y))
//Here we are cool with 1% accuracy.
#define FLOATS_APPROX_NEAR_TO(x,y) \
        BOOST_CHECK_MESSAGE( !isinf((x)) && !isinf((y))  && !isnan((x)) && !isnan((y)) \
                && fabs((x) - (y)) <= std::max<double>(fabs(x),fabs(y))/1e2 ,(x)<<" is not close to "<<(y))

BOOST_AUTO_TEST_CASE(check_sph_kern)
{
    //First check the kernel is working
    const double norm = 32./4/M_PI;
    FLOATS_NEAR_TO(sph_kernel(1),0);
    FLOATS_NEAR_TO(sph_kernel(0),norm);
    FLOATS_NEAR_TO(sph_kernel(0.5),0.25*norm);
    FLOATS_NEAR_TO(sph_kernel(0.25),0.71875*norm);
    FLOATS_NEAR_TO(sph_kernel(0.75),0.03125*norm);
}

BOOST_AUTO_TEST_CASE(check_sph_kern_frac)
{
   /* Integrals evaluated with Mathematica
    * If K(q) is the SPH kernel we need
    * int_zlow^zhigh K(q) dz
    *
    * Normalized such that 4 pi int_{q < 1} K(q) q^2 dq = 1
    * and q^2 = (b^2 + z^2)/h^2
    * The full range of the function is:
    * zlow = -1, zhigh = 1, coming to:
    *  2 int_0^1 K(z) dz = 3/(2 pi)
    *  zrange = sqrt(smooth*smooth-dr2)
    */
    //Correction from the normalisation I worked out in mathematica,
    //because the smoothing length definition changed.
    const double corr = 3/(4*M_PI);
    FLOATS_NEAR_TO(sph_cubic_kern_frac(-1,1,1,0,1), 3*2/M_PI);
    //Should be the same with a wide range.
    FLOATS_NEAR_TO(sph_cubic_kern_frac(-3,10,1,0,1), 3*2/M_PI);
    //Small variations
    FLOATS_APPROX_NEAR_TO(sph_cubic_kern_frac(-0.15,-0.1,1,0,1), corr*0.489167);
    //b!=0
    FLOATS_APPROX_NEAR_TO(sph_cubic_kern_frac(0.05,0.1,1,0.4,0.774597), corr*0.051004);
    FLOATS_APPROX_NEAR_TO(sph_cubic_kern_frac(0.15,0.16,1,0.8,0.447214), corr*0.000167423);
    FLOATS_APPROX_NEAR_TO(sph_cubic_kern_frac(-0.05,0.1,1,0.9,0.316228), corr*0.00040101);
    //Test wider range with b!=0
    FLOATS_APPROX_NEAR_TO(sph_cubic_kern_frac(0.3,1,1,0.3,0.83666), corr*0.177801);
    //Check outside of range
    FLOATS_NEAR_TO(sph_cubic_kern_frac(1.5,2,1,0,1), 0);


}


#define TNBINS 2000
BOOST_AUTO_TEST_CASE(check_compute_colden)
{
    //Give it the properties of Lyman alpha, a 10 Mpc box size, and a velfac from z=3.
    double velfac = 414.50523718485636/1e3 * 0.2 / 0.71;
    LineAbsorption test(1215.6701e-10,6.265e8,0.416400,1.00794,velfac, 10000, 1,SPH_CUBIC_SPLINE,1e-5);

    double colden[TNBINS] = {0};
    std::set<int> nonzero;
    /* Various checks on the column density*/
    //Add a particle with some reasonable properties that is exactly on the line of sight
    //halfway along the spectrum
    //Bins are 5kpc wide.
    test.add_colden_particle(colden, TNBINS, 0, 1,5002.5, 1);
    //Column density should be in one bin only
    //and total should be rho h int(-1,1)
    double total = 8*3./(4*M_PI);
    BOOST_CHECK_EQUAL(colden[999],0);
    FLOATS_NEAR_TO(colden[1000],total);
    BOOST_CHECK_EQUAL(colden[1001],0);
    nonzero.insert(1000);
    //Scales correctly with h:
    total /=4;
    //Change particle density with h
    test.add_colden_particle(colden, TNBINS, 0, 1/8.,1002.5, 2);
    BOOST_CHECK_EQUAL(colden[199],0);
    FLOATS_NEAR_TO(colden[200],total);
    BOOST_CHECK_EQUAL(colden[201],0);
    nonzero.insert(200);
    //Scales correctly with M:
    test.add_colden_particle(colden, TNBINS, 0, 10/8.,1012.5, 2);
    BOOST_CHECK_EQUAL(colden[201],0);
    FLOATS_NEAR_TO(colden[202],10*total);
    BOOST_CHECK_EQUAL(colden[203],0);
    nonzero.insert(202);
    //Split evenly over two bins
    test.add_colden_particle(colden, TNBINS, 0, 1/8.,1030, 2);
    BOOST_CHECK_EQUAL(colden[204],0);
    FLOATS_NEAR_TO(colden[205],total/2.);
    FLOATS_NEAR_TO(colden[206],total/2.);
    BOOST_CHECK_EQUAL(colden[207],0);
    nonzero.insert(205);
    nonzero.insert(206);
    //Correct periodic wrapping 
    test.add_colden_particle(colden, TNBINS, 0, 1/8.,10000, 2);
    BOOST_CHECK_EQUAL(colden[1998],0);
    FLOATS_NEAR_TO(colden[1999],total/2.);
    FLOATS_NEAR_TO(colden[0],total/2.);
    BOOST_CHECK_EQUAL(colden[1],0);
    nonzero.insert(1999);
    nonzero.insert(0);
    //Check non-zero offset from line
    test.add_colden_particle(colden, TNBINS, 0.7, 1,4852.5, 1);
    FLOATS_NEAR_TO(colden[969],0);
    FLOATS_APPROX_NEAR_TO(colden[970],0.0451531*3/(4*M_PI));
    BOOST_CHECK_EQUAL(colden[971],0);
    nonzero.insert(970);
    //Offset enough to not add anything
    test.add_colden_particle(colden, TNBINS, 1.0, 1,4802.5, 1);
    FLOATS_NEAR_TO(colden[960],0);

    //Check that all other values are zero still.
    for (int i=0; i<TNBINS;i++)
    {
        if(nonzero.count(i))
            continue;
        BOOST_CHECK_EQUAL(colden[i]+i,i);
    }

}

#define TLOS 25

//Macro to construct a sightline and store it in the next available place in the array
#define SIGHTLINE(a, x, y, z) do{ \
    cofm[3*nextlos]=(x); \
    cofm[3*nextlos+1]=(y); \
    cofm[3*nextlos+2]=(z); \
    axis[nextlos]=(a); \
    nextlos++;\
    } while(0)

BOOST_AUTO_TEST_CASE(check_index_table)
{
    int nextlos = 0;
    double cofm[TLOS*3];
    int axis[TLOS];
    //Construct a table of sightlines with various interesting and likely to break properties
    //Two sightlines with the same position.
    SIGHTLINE(1,4000,4000,4000);
    SIGHTLINE(1,4000,4000,4000);
    //Then one a little further away
    SIGHTLINE(1,4000,4020,4010);
    //Then one a little further in only one axis
    SIGHTLINE(1,4000,4000,4010);
    //One near the edge of the box
    SIGHTLINE(1,0,0.4,0.1);
    //Other box edge
    SIGHTLINE(1,0,10000-0.4,0.3);
    //One on its own
    SIGHTLINE(1,0,2000,1000);
    //Same but on different axes
    SIGHTLINE(3,1000,2000,500);
    SIGHTLINE(2,1000,2000,500);
    SIGHTLINE(1,1000,2000,500);
    //Fill out the table a bit
    SIGHTLINE(3,3000,5000,550);
    SIGHTLINE(3,8000,5500,9000);
    SIGHTLINE(1,6000,5500,9000);

    assert(nextlos < TLOS);
    //Construct the table
    IndexTable tab(cofm, axis, nextlos, 10000);

    //Now start testing.
    float pos[3] = {500,2000,1000};
    std::map<int, double> nearby = tab.get_near_lines(pos, 1);
    BOOST_CHECK_EQUAL(nearby.size(),1);
    BOOST_CHECK_EQUAL(nearby.begin()->first,6);
    BOOST_CHECK_EQUAL(nearby.begin()->second,0);
    //First coordinate not important
    float poss[3] = {5000,2000,1000};
    nearby = tab.get_near_lines(poss, 1);
    BOOST_CHECK_EQUAL(nearby.size(),1);
    BOOST_CHECK_EQUAL(nearby.begin()->first,6);
    BOOST_CHECK_EQUAL(nearby.begin()->second,0);
    //Slight offset, still within h
    float pos3[3] = {5010,2010,990};
    nearby = tab.get_near_lines(pos3, 20);
    BOOST_CHECK_EQUAL(nearby.size(),1);
    BOOST_CHECK_EQUAL(nearby.begin()->first,6);
    BOOST_CHECK_EQUAL(nearby.begin()->second,10*10+10*10.);
    //Slight offset, just outside h
    nearby = tab.get_near_lines(pos3, 10);
    BOOST_CHECK_EQUAL(nearby.size(),0);

    //Check duplicates are handled
    float pos2[3] = {4000,4000,4000};
    nearby = tab.get_near_lines(pos2, 1);
    BOOST_CHECK_EQUAL(nearby.size(),2);
    BOOST_CHECK_EQUAL(nearby.begin()->first,0);
    BOOST_CHECK_EQUAL(nearby.begin()->second,0);
    BOOST_CHECK_EQUAL((++nearby.begin())->first,1);
    BOOST_CHECK_EQUAL((++nearby.begin())->second,0);

    //Check duplicates are handled: wider
    nearby = tab.get_near_lines(pos2, 25);
    BOOST_CHECK_EQUAL(nearby.size(),4);
    BOOST_CHECK_EQUAL(nearby.at(0),0);
    BOOST_CHECK_EQUAL(nearby.at(1),0);
    BOOST_CHECK_EQUAL(nearby.at(2),20*20+10*10.);
    BOOST_CHECK_EQUAL(nearby.at(3),10*10.);

    //Check periodic wrapping is working
    float pos5[3] = {1000,9999.9,9999.9};
    nearby = tab.get_near_lines(pos5, 0.6);
    BOOST_CHECK_EQUAL(nearby.size(),2);
    FLOATS_APPROX_NEAR_TO(nearby.at(4),0.5*0.5+0.2*0.2);
    FLOATS_APPROX_NEAR_TO(nearby.at(5),0.3*0.3+0.4*0.4);

    //Check multiple axes
    float pos4[3] = {1000.5,2000,501};
    nearby = tab.get_near_lines(pos4, 1.5);
    BOOST_CHECK_EQUAL(nearby.size(),3);
    FLOATS_APPROX_NEAR_TO(nearby.at(7),0.25);
    FLOATS_APPROX_NEAR_TO(nearby.at(8),0.25+1);
    FLOATS_APPROX_NEAR_TO(nearby.at(9),1);

    //Now test get_near_particles
    float poses[3*9] = { 500,2000,1000, 5000,2000.0,1000.0, 5010.0,2010,990,
                        4000,4000,4000, 1000,9999.9,9999.9, 1000.5,2000,501,
                        7500,7500,7500, 4008,4008.0,4008.0, 2000.0,9999,9999.8};
    float hh[9] = {1,1,20,25,0.6,1.5,7,10,0.8};
    std::valarray< std::map<int, double> > nearby_array = tab.get_near_particles(poses, hh, 9);
    BOOST_CHECK_EQUAL(nearby_array.size(), nextlos);
    //Did we pick up the right number of particles in all cases?
    BOOST_CHECK_EQUAL(nearby_array[0].size(),1);
    BOOST_CHECK_EQUAL(nearby_array[1].size(),1);
    BOOST_CHECK_EQUAL(nearby_array[2].size(),1);
    BOOST_CHECK_EQUAL(nearby_array[3].size(),2);
    BOOST_CHECK_EQUAL(nearby_array[4].size(),1);
    BOOST_CHECK_EQUAL(nearby_array[5].size(),2);
    BOOST_CHECK_EQUAL(nearby_array[6].size(),3);
    BOOST_CHECK_EQUAL(nearby_array[7].size(),1);
    BOOST_CHECK_EQUAL(nearby_array[8].size(),1);
    BOOST_CHECK_EQUAL(nearby_array[9].size(),1);
    BOOST_CHECK_EQUAL(nearby_array[10].size(),0);
    BOOST_CHECK_EQUAL(nearby_array[11].size(),0);
    BOOST_CHECK_EQUAL(nearby_array[12].size(),0);
    //Check a few values
    BOOST_CHECK_EQUAL(nearby_array[0].begin()->first,3);
    BOOST_CHECK(nearby_array[3].find(3) != nearby_array[3].end());
    BOOST_CHECK(nearby_array[3].find(3) != nearby_array[3].end());
    BOOST_CHECK(nearby_array[5].find(4) != nearby_array[5].end());
    BOOST_CHECK(nearby_array[5].find(8) != nearby_array[5].end());
    std::map<int, double>::iterator it = nearby_array[6].begin();
    BOOST_CHECK_EQUAL(it->first,0);
    BOOST_CHECK_EQUAL((++it)->first,1);
    BOOST_CHECK_EQUAL((++it)->first,2);
    BOOST_CHECK_EQUAL(nearby_array[8].begin()->first,5);
}

BOOST_AUTO_TEST_CASE(check_profile)
{
    //Check that the Voigt profile is as we expect.
    //Integrals evaluated exactly with mathematica.
    FLOATS_NEAR_TO(profile(0,0),1);
    FLOATS_NEAR_TO(profile(0,0.1),0.896457);
    //Close to center
    FLOATS_NEAR_TO(profile(0.1,1e-6),0.990048);
    FLOATS_NEAR_TO(profile(0.1,1e-4),0.989939);
    //Far from center
    FLOATS_NEAR_TO(profile(15,1e-6),2.52441e-9);
    FLOATS_NEAR_TO(profile(15,1e-4),2.52441e-7);
    FLOATS_NEAR_TO(profile(20,1e-7),1.4158e-10);
    //Intermediate
    FLOATS_NEAR_TO(profile(1,1e-4),0.367888);
    FLOATS_NEAR_TO(profile(1.5,1e-6),0.1054);
    FLOATS_NEAR_TO(profile(2,1e-7),0.0183157);
    FLOATS_NEAR_TO(profile(1,1e-3),0.367965);
}

BOOST_AUTO_TEST_CASE(check_single_absorber)
{
    //Integrals evaluated numerically as usual with mathematica
    //H at 20000 K.
    double bb = 0.128557*sqrt(2e4/1);
    //Correction from the normalisation I worked out in mathematica,
    //because the smoothing length definition changed.
    const double corr = 3/(4*M_PI);
    //First test not offset
    SingleAbsorber sing(bb,0,10,1.e-4, SPH_CUBIC_SPLINE);
    //Absorption at the origin
    FLOATS_NEAR_TO(sing.tau_kern_outer(0,0),corr*78.0409);
    FLOATS_NEAR_TO(sing.tau_kern_outer(5,5),corr*72.6216);
    FLOATS_NEAR_TO(sing.tau_kern_outer(10,10),corr*58.5185);
    FLOATS_NEAR_TO(sing.tau_kern_outer(20,20),corr*24.6696);
    //Try different smoothing length
    SingleAbsorber sing2(bb,0,2,1.e-4, SPH_CUBIC_SPLINE);
    FLOATS_NEAR_TO(sing2.tau_kern_outer(5,5),corr*14.8203);
    //Offset from the sightline
    SingleAbsorber sing3(bb,25,10,1.e-4, SPH_CUBIC_SPLINE);
    FLOATS_APPROX_NEAR_TO(sing3.tau_kern_outer(0,0),corr*18.218);
    FLOATS_APPROX_NEAR_TO(sing3.tau_kern_outer(5,5),corr*16.9436);
    //Check integrating over a range in v
    //Carbon
    bb = 0.128557*sqrt(2e4/16);
    SingleAbsorber sing4(bb,0,5,1.e-6, SPH_CUBIC_SPLINE);
    FLOATS_APPROX_NEAR_TO(sing4.tau_kern_outer(0,10),corr*16.0403);
    FLOATS_APPROX_NEAR_TO(sing4.tau_kern_outer(-5,5),corr*27.1978);
}

#define  BOLTZMANN  1.3806504e-16  /*  ergs K-1 or cm2 g s-2 K-1 */
#define  PROTONMASS  1.67262178e-24 /* 1 a.m.u */
//Little macros to help in the below
#define vbin(bin, pos) (10000*velfac/TNBINS*(bin) - velfac*(pos))
#define SA(bin, pos, dens) (amp*(dens)/bb*sing.tau_kern_outer(vbin((bin),(pos)), vbin((bin+1),(pos)))/velfac)
BOOST_AUTO_TEST_CASE(check_add_tau)
{
    //Give it the properties of Lyman alpha, a 10 Mpc box size, and a velfac from z=3.
    double velfac = 414.50523718485636/1e3 * 0.2 / 0.71;
    LineAbsorption test(1215.6701e-10,6.265e8,0.416400,1.00794,velfac, 10000, 0.25, SPH_CUBIC_SPLINE,1e-5);
    //Since SingleAbsorber is tested above, use it in this test
    double temp = 2e4;
    double bb  = sqrt(2.0*BOLTZMANN/(PROTONMASS))/1e5*sqrt(temp/1.00794);
    double smooth = 3;
    double amp = 7.57973e-15;
    double voigt_fac = 1215.6701e-10*6.265e8/(4.*M_PI)/1e5;
    SingleAbsorber sing(bb,0,velfac*smooth,voigt_fac/bb, SPH_CUBIC_SPLINE);
    //Conversion factor from cm to kpc/h at z=3
    double rscale = 3.085678e21*0.25/0.7;
    double tau[TNBINS] = {0};
    /* Various checks on the column density*/
    //Bins are 5kpc (1.5km/s) wide.
    test.add_tau_particle(tau, TNBINS, 0, 1e-3*rscale,5002.5, 0,temp, smooth);
    //Compare resulting tau to explicit calls to SingleAbsorber
    FLOATS_NEAR_TO(tau[999],SA(999,5002.5, 1e-3*rscale));
    FLOATS_NEAR_TO(tau[1000],SA(1000, 5002.5,1e-3*rscale));
    //Symmetric function
    FLOATS_NEAR_TO(tau[1001],tau[999]);
    //Check tail
    FLOATS_NEAR_TO(tau[400],SA(400,5002.5, 1e-3*rscale));
    FLOATS_NEAR_TO(tau[401],SA(401,5002.5, 1e-3*rscale));
    FLOATS_NEAR_TO(tau[402],SA(402,5002.5, 1e-3*rscale));
    //Reset tau
    memset(tau, 0, TNBINS*sizeof(double));
    //Check peculiar vel dependence: this should just induce an offset
    //3 bins
    double pecvel = 10000/TNBINS*3*velfac/sqrt(0.25);
    test.add_tau_particle(tau, TNBINS, 0, 1e-3*rscale,5002.5,pecvel ,temp, smooth);
    FLOATS_NEAR_TO(tau[1002],SA(999,5002.5, 1e-3*rscale));
    FLOATS_NEAR_TO(tau[1003],SA(1000, 5002.5,1e-3*rscale));
    //Symmetric function
    FLOATS_NEAR_TO(tau[1004],tau[1002]);
    memset(tau, 0, TNBINS*sizeof(double));
    //Offset enough to not add anything
    test.add_tau_particle(tau, TNBINS, 1, 1,5002.5, 0,10000, 1);
    FLOATS_NEAR_TO(tau[1000],0);
    FLOATS_NEAR_TO(tau[1501],0);
}

//Check that if we compute column density and optical depth with the same input, we get a consistent answer.
BOOST_AUTO_TEST_CASE(check_tau_colden_consistency)
{
    //Give it the properties of Lyman alpha, a 10 Mpc box size, and a velfac from z=3.
    double velfac = 414.50523718485636/1e3 * 0.2 / 0.71;
    LineAbsorption test(1215.6701e-10,6.265e8,0.416400,1.00794,velfac, 10000, 0.25,SPH_CUBIC_SPLINE,1e-5);
    //Since SingleAbsorber is tested above, use it in this test
    double temp = 2e4;
    double smooth = 3;
    double amp = 7.57973e-15;
    //Conversion factor from cm to kpc/h at z=3
    double rscale = 3.085678e21*0.25/0.7;
    double tau[TNBINS] = {0};
    double colden[TNBINS] = {0};
    //Bins are 5kpc (1.5km/s) wide: add particle to both tau and colden.
    //Check that the ratio between optical depth and column density is correct.
    test.add_tau_particle(tau, TNBINS, 0, 0.1*rscale,5002.5, 0,temp, smooth);
    test.add_colden_particle(colden, TNBINS, 0, 0.1*rscale,5002.5, smooth);
    double scolden = std::accumulate(colden, colden+TNBINS,0.);
    double stau =  std::accumulate(tau, tau+TNBINS,0.);
    FLOATS_NEAR_TO(stau/scolden, amp/velfac/2.81809);
}
