#ifndef FAST_VOIGT_H
#define FAST_VOIGT_H

#include <cmath>
#include "Faddeeva.h"

/* A fast evaluation of the Voigt function H(a,u) = Re[w(u + i a)], which is
 * where almost all of the time in compute_tau goes.
 *
 * The damping parameter a = gamma lambda / (4 pi b) is small: 5e-4 for
 * Lyman alpha at 10^4 K, and below 3e-3 for any line and temperature we
 * encounter. In that regime the expansion in a (Harris 1948) is essentially
 * exact:
 *
 *      H(a,u) = exp(-u^2) + a H1(u) + O(a^2)
 *      H1(u)  = (2/sqrt(pi)) (2 u Daw(u) - 1) = 2 u Im[w(u)] - 2/sqrt(pi)
 *
 * The neglected term is a^2 (1 - 2u^2) exp(-u^2), so for a = 5e-4 the error
 * is 2.5e-7 absolute and 6e-6 relative. Above VOIGT_AMAX we call Faddeeva
 * rather than extrapolate the expansion into a regime it is not good for.
 *
 * H1 is smooth, O(1), and multiplied by a small number, so tabulating it and
 * interpolating linearly leaves an error well below the a^2 truncation.
 * The Gaussian is tabulated on the same grid: exp() is otherwise the most
 * expensive thing left once Faddeeva is gone.
 *
 * Past VOIGT_UMAX the Gaussian has underflowed and only the asymptotic series
 * for H1 is needed. That is where the damping wings are, which is where most
 * of the bins are for a strong absorber, so those bins cost a few flops. */

#define VOIGT_UMAX 6.0
#define VOIGT_NTAB 1024
/* Above this a the expansion is not accurate enough and we call Faddeeva. */
#define VOIGT_AMAX 0.001

struct VoigtTable
{
    /*The Gaussian, the first order damping correction, and the difference
     * between neighbouring entries of the latter (so interpolating is a
     * multiply-add rather than two loads and a subtraction).*/
    double gau[VOIGT_NTAB+2];
    double h1[VOIGT_NTAB+2];
    double dh1[VOIGT_NTAB+2];
    VoigtTable()
    {
        for(int i = 0; i < VOIGT_NTAB+2; ++i) {
            const double u = i*(VOIGT_UMAX/VOIGT_NTAB);
            gau[i] = exp(-u*u);
            h1[i] = 2*u*Faddeeva::w_im(u) - 2/sqrt(M_PI);
        }
        for(int i = 0; i < VOIGT_NTAB+1; ++i)
            dh1[i] = h1[i+1] - h1[i];
        dh1[VOIGT_NTAB+1] = 0;
    }
};

/*Built once, at startup: defined in absorption.cpp.*/
extern const VoigtTable voigt_table;

/* H(a,u) to first order in a. Only call this for a < VOIGT_AMAX;
 * profile() in singleabs.h does the check. */
inline double fast_voigt(const double uu, const double aa)
{
    const double u = fabs(uu);
    if (u >= VOIGT_UMAX) {
        /*Asymptotic series for H1: the Gaussian is < 3e-16 here.*/
        const double u2 = 1./(u*u);
        return aa*M_2_SQRTPI/2*u2*(1 + u2*(1.5 + u2*(3.75 + u2*(13.125 + 59.0625*u2))));
    }
    const double x = u*(VOIGT_NTAB/VOIGT_UMAX);
    const int i = (int) x;
    const double f = x - i;
    /*exp(-u^2) = exp(-u_i^2) exp(-(2 u_i d + d^2)) where d is the distance to
     * the tabulated point below. The second factor is within 0.08 of one, so
     * a short Taylor series for it is good to 1e-8.*/
    const double d = f*(VOIGT_UMAX/VOIGT_NTAB);
    const double t = -d*(2*i*(VOIGT_UMAX/VOIGT_NTAB) + d);
    const double ex = voigt_table.gau[i]*(1+t*(1+t*(0.5+t*(1./6+t*(1./24)))));
    return ex + aa*(voigt_table.h1[i] + f*voigt_table.dh1[i]);
}

#endif
