#ifndef ABSORPTION_H
#define ABSORPTION_H

#define  PROTONMASS  1.67262178e-24 /* 1 a.m.u */

/* Class to compute the absorption from a single particle
 * onto a spectrum
 * */
class LineAbsorption
{
    public:
        /*Dimensions:
         * velocities in km/s (physical).
         * distances in kpc/h (comoving)
         * velfac: factor to convert from distance to velocity units.
         * Should be  atime /h100 * Hz/1e3 (Hz in km/s/Mpc)
         * boxsize: Size of the box in comoving kpc/h
         * atime: Scale factor
         * */
        LineAbsorption(const double lambda, const double gamma, const double fosc, const double amumass, const double velfac_i, const double boxsize, const double atime_i);

        /* Add the absorption from a particle to the spectrum in the array
         * tau, and the density from the particle to the array colden
         * The slightly C-style interface is so we can easily use the data in python.
         *
         * Output:
         * tau: array specifying the optical depth of the spectrum.
         * If this is NULL, just compute the column density.
         * colden: array specifying the column density of the spectrum. (1e10 M_sun/h / (comoving kpc/h)^2)
         * nbins: Size of above arrays
         *
         * Input:
         * dr2: transverse distance to spectra from particle (comoving kpc/h)
         * dens: value of density field for absorbing species (1e10 M_sun / h / (comoving kpc/h)^3 or amu/cm^3)
         * ppos: particle distance from box edge parallel to spectrum (comoving kpc/h)
         * pvel: particle velocity parallel to spectrum (physical km/s)
         * temp: particle temperature (K)
         * smooth: particle smoothing length (comoving kpc/h)
         */
        void add_particle(double * tau, double * colden, const int nbins, const double dr2, const float dens, const float ppos, const float pvel, const float temp, const float smooth);

        /* Compute the absorption in a single bin, using
         * either straight Gaussian or a Voigt profile.
         * Note! No unit conversion is done on the column density.
         * Arguments:
         * colden: column density of absorber in amu per m^2.
         * vdiff: the relative velocities between absorper and bin.
         * b_H1: b parameter of absorber. = bfac * sqrt(temp), temperature of absorber in K
         */
        inline double tau_single(const double colden, const double vdiff, const double b_H1);

    private:
        /* Absorption cross-sections m^2 */
        const double sigma_a;
        /* Constant factor to turn sqrt(temperature) into velocity*/
        const double bfac;
        /* Factor to turn b into a dimensionless Voigt profile broadening factor,
         * giving the balance between doppler and thermal broadening. */
        const double voigt_fac;
        const double velfac, vbox, atime;
};

/* Compute temperature (in K) from internal energy.
 * uu: internal energy in Gadget units
 * ne: electron abundance
 * xh: hydrogen mass fraction (0.76)
 * Factor to convert U (J/kg) to T (K) : U = N k T / (γ - 1)
 * T = U (γ-1) μ m_P / k_B
 * where k_B is the Boltzmann constant
 * γ is 5/3, the perfect gas constant
 * m_P is the proton mass
 * μ is 1 / (mean no. molecules per unit atomic weight)
 */
double compute_temp(const double uu, const double ne, const double xh);

#endif
