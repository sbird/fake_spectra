#ifndef ABSORPTION_H
#define ABSORPTION_H

#define  PROTONMASS  1.67262178e-24 /* 1 a.m.u */

/*Type for the kernel function pointer*/
typedef double (*const kern_frac_func)(double, double, const double, const double, const double);

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
         * lambda: wavelength in cm
         * gamma: transition rate in 1/s
         * fosc: oscillation fraction
         * */
        LineAbsorption(const double lambda, const double gamma, const double fosc, const double amumass, const double velfac_i, const double boxsize, const double atime_i, const int kernel_i, const double tautail);

        /* Add the absorption from a particle to the spectrum in the array
         * tau, or the density from the particle to the array colden,
         * depending on function called.
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
         * dens: value of density field for absorbing species (1e10 M_sun / h / (comoving kpc/h)^3 from C extractor 
         * or amu/cm^2 kpc/h from the python module)
         * ppos: particle distance from box edge parallel to spectrum (comoving kpc/h)
         * pvel: particle velocity parallel to spectrum (physical km/s) (not needed for colden)
         * temp: particle temperature (K) (not needed for colden)
         * smooth: particle smoothing length (comoving kpc/h)
         */
        void add_colden_particle(double * colden, const int nbins, const double dr2, const float dens, const float ppos, const float smooth);
        void add_tau_particle(double * tau, const int nbins, const double dr2, const float dens, const float ppos, const float pvel, const float temp, const float smooth);

    private:
        /* Threshold of tau below which we stop computing profiles.
         * Note that because this is for each particle, it should be fairly small.*/
        const double tautail;
        /* Absorption cross-sections cm^2 */
        const double sigma_a;
        /* Constant factor to turn sqrt(temperature) in K into velocity in km/s*/
        const double bfac;
        /* Factor to turn b into a dimensionless Voigt profile broadening factor,
         * giving the balance between doppler and thermal broadening in km/s. */
        const double voigt_fac;
        const double velfac, vbox, atime;
        const kern_frac_func kern_frac;
        const int kernel;
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
