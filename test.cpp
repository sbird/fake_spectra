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


BOOST_AUTO_TEST_CASE(check_index_table)
{
}
