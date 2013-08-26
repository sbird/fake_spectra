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
#include "index_table.h"
#include <cstdio>
#include <math.h>
#include <set>
#include <boost/test/unit_test.hpp>
#include <boost/test/test_tools.hpp>

#define FLOATS_NEAR_TO(x,y) \
        BOOST_CHECK_MESSAGE( !isinf((x)) && !isinf((y))  && !isnan((x)) && !isnan((y)) \
                && fabs((x) - (y)) <= std::max<float>(fabs(x),fabs(y))/1e5 ,(x)<<" is not close to "<<(y))
//Here we are cool with 1% accuracy.
#define FLOATS_APPROX_NEAR_TO(x,y) \
        BOOST_CHECK_MESSAGE( !isinf((x)) && !isinf((y))  && !isnan((x)) && !isnan((y)) \
                && fabs((x) - (y)) <= std::max<float>(fabs(x),fabs(y))/1e2 ,(x)<<" is not close to "<<(y))

double sph_kern_frac(double zlow, double zhigh, double bb2);

BOOST_AUTO_TEST_CASE(check_sph_kern)
{
    /* Integrals evaluated exactly with Mathematica*/
    /* If K(q) is the SPH kernel normalized such that 4 pi int_{q < 1} K(q) q^2 dq = 1
    * and q^2 = b^2 + z^2, then sph_kern_frac computes:
    * int_zlow^zhigh K(q) dz.
    * The full range of the function is:
    * zlow = -1, zhigh = 1, coming to:
    *  2 int_0^1 K(z) dz = 6 / pi
    */
    FLOATS_NEAR_TO(sph_kern_frac(-1,1,0), 6./M_PI);
    //Should be the same with a wide range.
    FLOATS_NEAR_TO(sph_kern_frac(-3,10,0), 6./M_PI);
    //Small variations
    FLOATS_APPROX_NEAR_TO(sph_kern_frac(-0.15,-0.1,0), 0.11678);
    //b!=0
    FLOATS_APPROX_NEAR_TO(sph_kern_frac(0.05,0.1,0.4), 0.0121763);
    FLOATS_APPROX_NEAR_TO(sph_kern_frac(0.15,0.16,0.8), 3.99694e-5);
    FLOATS_APPROX_NEAR_TO(sph_kern_frac(-0.05,0.1,0.9), 9.57342e-5);
    //Test wider range with b!=0
    FLOATS_APPROX_NEAR_TO(sph_kern_frac(0.3,1,0.3), 0.0424468);
}


#define TNBINS 2000
BOOST_AUTO_TEST_CASE(check_compute_colden)
{
    //Give it the properties of Lyman alpha, a 10 Mpc box size, and a velfac from z=3.
    double velfac = 414.50523718485636/1e3 * 0.2 * 0.71;
    ComputeLineAbsorption test(1215.6701e-10,6.265e8,0.416400,1.00794,velfac, 10000);

    double* tau = NULL;
    double colden[TNBINS] = {0};
    std::set<int> nonzero;
    /* Various checks on the column density*/
    //Add a particle with some reasonable properties that is exactly on the line of sight
    //halfway along the spectrum, not moving and has a temperature of 10^4 K.
    //Bins are 5kpc wide.
    test.add_particle(tau, colden, TNBINS, 0, 1,5002.5,0, 10000, 1);
    //Column density should be in one bin only
    //and total should be dz*M/(h^3) int(-1,1)
    BOOST_CHECK_EQUAL(colden[999],0);
    FLOATS_NEAR_TO(colden[1000],10000./TNBINS*6/M_PI);
    BOOST_CHECK_EQUAL(colden[1001],0);
    nonzero.insert(1000);
    //Scales correctly with h:
    test.add_particle(tau, colden, TNBINS, 0, 1,1002.5,0, 10000, 2);
    BOOST_CHECK_EQUAL(colden[199],0);
    FLOATS_NEAR_TO(colden[200],10000./TNBINS*6/M_PI/pow(2,3));
    BOOST_CHECK_EQUAL(colden[201],0);
    nonzero.insert(200);
    //Scales correctly with M:
    test.add_particle(tau, colden, TNBINS, 0, 10,1012.5,0, 10000, 2);
    BOOST_CHECK_EQUAL(colden[201],0);
    FLOATS_NEAR_TO(colden[202],10*10000./TNBINS*6/M_PI/pow(2,3));
    BOOST_CHECK_EQUAL(colden[203],0);
    nonzero.insert(202);
    //Split evenly over two bins
    test.add_particle(tau, colden, TNBINS, 0, 1,1030,0, 10000, 2);
    BOOST_CHECK_EQUAL(colden[204],0);
    FLOATS_NEAR_TO(colden[205],10000./TNBINS*6/M_PI/pow(2,3)/2.);
    FLOATS_NEAR_TO(colden[206],10000./TNBINS*6/M_PI/pow(2,3)/2.);
    BOOST_CHECK_EQUAL(colden[207],0);
    nonzero.insert(205);
    nonzero.insert(206);
    //Correct periodic wrapping 
    test.add_particle(tau, colden, TNBINS, 0, 1,10000,0, 10000, 2);
    BOOST_CHECK_EQUAL(colden[1998],0);
    FLOATS_NEAR_TO(colden[1999],10000./TNBINS*6/M_PI/pow(2,3)/2.);
    FLOATS_NEAR_TO(colden[0],10000./TNBINS*6/M_PI/pow(2,3)/2.);
    BOOST_CHECK_EQUAL(colden[1],0);
    nonzero.insert(1999);
    nonzero.insert(0);
    //Correct shifting due to peculiar velocity:
    //velfac is 0.0589, so this is a shift by 50kpc = 10 bins.
    test.add_particle(tau, colden, TNBINS, 0, 1,5000,2.9429871840124799, 10000, 2);
    BOOST_CHECK_EQUAL(colden[1008],0);
    FLOATS_NEAR_TO(colden[1009],10000./TNBINS*6/M_PI/pow(2,3)/2.);
    FLOATS_NEAR_TO(colden[1010],10000./TNBINS*6/M_PI/pow(2,3)/2.);
    BOOST_CHECK_EQUAL(colden[1011],0);
    nonzero.insert(1009);
    nonzero.insert(1010);
    //Correct shifting due to negative peculiar velocity:
    //velfac is 0.0589, so this is a shift by -50kpc = -10 bins.
    test.add_particle(tau, colden, TNBINS, 0, 1,5000,-2.9429871840124799, 10000, 2);
    BOOST_CHECK_EQUAL(colden[988],0);
    FLOATS_NEAR_TO(colden[989],10000./TNBINS*6/M_PI/pow(2,3)/2.);
    FLOATS_NEAR_TO(colden[990],10000./TNBINS*6/M_PI/pow(2,3)/2.);
    BOOST_CHECK_EQUAL(colden[991],0);
    nonzero.insert(989);
    nonzero.insert(990);
    //Correct periodic wrapping due to velocity shifting
    //velfac is 0.0589, so this is a shift by -20Mpc - 100kpc = -20 bins.
    test.add_particle(tau, colden, TNBINS, 0, 1,5000,-1183.0808479730169, 10000, 2);
    BOOST_CHECK_EQUAL(colden[978],0);
    FLOATS_APPROX_NEAR_TO(colden[979],10000./TNBINS*6/M_PI/pow(2,3)/2.);
    FLOATS_APPROX_NEAR_TO(colden[980],10000./TNBINS*6/M_PI/pow(2,3)/2.);
    BOOST_CHECK_EQUAL(colden[981],0);
    nonzero.insert(979);
    nonzero.insert(980);
    //Check non-zero offset from line
    test.add_particle(tau, colden, TNBINS, 0.7, 1,4852.5,0, 10000, 1);
    FLOATS_APPROX_NEAR_TO(colden[969],0);
    FLOATS_APPROX_NEAR_TO(colden[970],10000./TNBINS*0.0107795);
    BOOST_CHECK_EQUAL(colden[971],0);
    nonzero.insert(970);
    //Offset enough to not add anything
    test.add_particle(tau, colden, TNBINS, 1.0, 1,4802.5,0, 10000, 1);
    FLOATS_APPROX_NEAR_TO(colden[960],0);

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
    std::valarray< std::vector<std::pair <int, double> > > nearby_array = tab.get_near_particles(poses, hh, 9);
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
    BOOST_CHECK_EQUAL(nearby_array[3][0].first,3);
    BOOST_CHECK_EQUAL(nearby_array[3][1].first,7);
    BOOST_CHECK_EQUAL(nearby_array[5][0].first,4);
    BOOST_CHECK_EQUAL(nearby_array[5][1].first,8);
    std::vector<std::pair <int, double> >::iterator it = nearby_array[6].begin();
    BOOST_CHECK_EQUAL(it->first,0);
    BOOST_CHECK_EQUAL((++it)->first,2);
    BOOST_CHECK_EQUAL((++it)->first,1);
    BOOST_CHECK_EQUAL(nearby_array[8][0].first,5);
}
