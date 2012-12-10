"""Python module for generating fake spectra from an N-body catalogue."""

import numpy as np
import hsml
from _spectra_priv import _SPH_Interpolate

def SPH_Interpolate(data, los_table, nbins, box):
    """Interpolate particles to lines of sight, calculating density, temperature and velocity
    of various species (TODO: only hydrogen now) along the line of sight.

    This is a wrapper which calls the C function.
    Arguments:
    	data - HDF5 dataset from snapshot. Use f["PartType0"]
	    los_table - table of los positions. should have member arrays x, y, z and axis.
	    nbins - number of bins in each spectrum
	    box - box size

    Returns:
        rho_HI
        vel_HI
        temp_HI
        all as arrays along the line of sight
    """
    pos = np.array(data["Coordinates"],dtype=np.float32)
    vel = np.array(data["Velocities"],dtype=np.float32)
    mass = np.array(data["Masses"],dtype=np.float32)
    u = np.array(data["InternalEnergy"],dtype=np.float32)
    nh0 = np.array(data["NeutralHydrogenAbundance"],dtype=np.float32)
    ne = np.array(data["ElectronAbundance"],dtype=np.float32)
    try:
        metals = np.array(data["GFM_Metals"],dtype=np.float32)[:,2:]
    except IOError:
        metals = np.array()
    hh = np.array(hsml.get_smooth_length(data),dtype=np.float32)
    xx=np.array(los_table.xx, dtype=np.float32)
    yy=np.array(los_table.yy, dtype=np.float32)
    zz=np.array(los_table.zz, dtype=np.float32)
    axis=np.array(los_table.axis, dtype=np.int32)
    return  _SPH_Interpolate(nbins, box, pos, vel, mass, u, nh0, ne, metals, hh, axis, xx, yy, zz)


#Speed of light
C = 2.99e8
#Boltzmann constant
KBOLTZ =1.3806504e-23

class line:
    sigma_X
    gamma_X
    lambda_X
    m_X

"""Note in Arepo we have GFM_Metals and GFM_Metallicity.

GFM_Metallicity is the total mass in species not H or He
per unit gas mass (and is used for cooling).

GFM_Metals is a 9-component array of species:
H, He, C, N, O, Ne, Mg, Si, Fe

Because these are not all the species, GFM_Metals will not sum to 1
and sum(GFM_Metals[2:])  < GFM_Metallicity

However, it should be true that
1- sum(GFM_Metals[:]) + GFM_Metals[2:] = GFM_Metallicity
"""

def compute_absorption(xbins, rho, vel, temp, line):
    """Computes the absorption spectrum (tau (u) ) from a binned set of interpolated densities, velocities and temperatures.
    xbins are the positions of each bin along the sightline. A good default is
    xbins = np.range(0,NBINS)*Box/NBINS
    Adapted from ComputeAbsorption in flux_extractor.
    sigma_X is the cross-section for this transition.
    Optical depth is given by:
        \\tau (u) = \sigma_X c / H(z) \int_\infty^\infty n_x(x) V( u - x - v_pec, b(x) ) dx
        where V is the Voigt profile, b(x)^2 = 2k_B T /m_x c^2 is the velocity dispersion.
        and v_pec is the peculiar velocity.
        """
#  const double H0 = 1.0e5/MPC; /* 100kms^-1Mpc^-1 in SI */
#    /* Critical matter/energy density at z = 0.0 */
#  const double rhoc = 3.0 * (H0*h100)*(H0*h100) / (8.0 * M_PI * GRAVITY); /* kgm^-3 */
#  /* Mean hydrogen mass density of the Universe */
#  const double critH = (rhoc * omegab * XH) / (atime*atime*atime); /* kgm^-3*/
#  /* Conversion factors from internal units */
#  const double rscale = (KPC*atime)/h100;   /* convert length to m */
#  const double vscale = sqrt(atime);        /* convert velocity to kms^-1 */
#  const double mscale = (1.0e10*SOLAR_MASS)/h100; /* convert mass to kg */
#  const double escale = 1.0e6;           /* convert energy/unit mass to J kg^-1 */
#  const double tscale = ((GAMMA-1.0) * HMASS * PROTONMASS * escale ) / BOLTZMANN; /* convert (with mu) T to K */
#  /*const double hscale = rscale * 0.5;*/ /* Note the factor of 0.5 for this kernel definition */
#  /*    Calculate the length scales to be used in the box */
#  const double dzgrid   = (box100) / (double)NBINS; /* bin size (kpc) */
#  const double vmax = box100 * Hz * rscale/ MPC; /* box size (kms^-1) */
#  const double vmax2 = vmax/2.0; /* kms^-1 */
#  const double dvbin = vmax / (double)NBINS; /* bin size (kms^-1) */

#   const double sigma_Lya_H1  = sqrt(3.0*M_PI*SIGMA_T/8.0) * LAMBDA_LYA_H1  * FOSC_LYA;
  #Optical depth prefactor
  tau_prefac = rscale*line.sigma_X*C*dzgrid/sqrt(math.pi)
  # If there are no particles in this bin, rhoker will be zero.
  # In this case, we set temp and veloc arbitrarily to one,
  # to avoid nans propagating. Zero rho will imply zero absorption
  # anyway.
  #ind = np.where(rho > 0)
  #veloc[ind]  *= vscale/rho[ind] # HI weighted km s^-1
  #temp[ind]   *= tscale/rho[ind # HI weighted K
  #rho *= mscale/rscale**3   #Put rho in m units
  #ind = np.where(rho <=0)
  #veloc[ind]=1
  #temp[ind]=1
  #Compute the spectra
  tau = np.zeros(NBINS)
  #Total velocity of each bin relative to start of box: H x + v
  u_H  = Hubble(z) /MPC *xbins + vel
  vmax  = Hubble(z) /MPC *xbins[-1]
  for i in xrange(0, NBINS):
      #Velocity of this cloud relative to this bin, wrapped periodically
      vdiff  = fabs(xbins[i] - u_H1)
      ind = np.where(vdiff > vmax)
      vdiff[ind] = vmax - vdiff
      bb = vel_disp(temp, line.m_x)
      #x for the Voigt profile is (c/b) (lambda/lambda_i-1), so this is x/c
      x = (vdiff/bb
      #a for the Voigt profile
      aa = line.gamma_X*line.lambda_X/(4.0*math.pi*bb)
      prof = voigt(aa, x)
      tau[i] = A_H1 * np.sum(rho*profile/bb/m_x)
  return tau

def vel_disp(temp, m_x):
    """Velocity dispersion (b param) as a function of temperature and particle mass"""
    return np.sqrt(2*KBOLTZ*temp/m_x/C**2)

def voigt(a, x):
    """Voigt profile (Tepper-Garcia, 2006: astro-ph/0602124)"""
    x2 = x*x
    T1 = exp(-x2)
    T2 = 1.5/x2
    if x2 < 1.0e-10:
        profile = T1
    else:
        profile = T1 - a/sqrt(math.pi)/x2*(T1*T1*(4.0*x2**2 + 7.0*x2 + 4.0 + T2) - T2 -1.0)
    return profile

#Stolen from yt
def voigt(a,u):
    """
    NAME:
        VOIGT
    PURPOSE:
        Implementation of Voigt function
    CATEGORY:
            Math
    CALLING SEQUENCE:
            voigt=Voigt(a,u)
    INPUTS:
            A = Voigt "A" parameter.
            U = Frequency in units of the Doppler frequency.

            The line profile "Phi(v)", the doppler width
            "Delv", the voigt parameter "a", and the frequency "u"
            are given by:

            Phi(v) =  Voigt(a,u)/[ Delv * sqrt(pi) ]
            Delv   =  Vo/c * sqrt[ 2kT/m ]
            u      =  V - Vo / Delv
            a      =  GAMMA / [ Delv * 4pi ]
            Gamma  =  Gu + Gl + 2*Vcol
            "Gu" and "Gl" are the widths of the upper and lower states
            "Vcol" is the collisions per unit time
            "Vo" is the line center frequency

    OUTPUTS:
            An array of the same type as u
    RESTRICTIONS:
            U must be an array, a should not be. Also this procedure is only valid
            for the region a<1.0, u<4.0 or a<1.8(u+1), u>4, which should be most
            astrophysical conditions (see the article below for further comments
    PROCEDURE:
            Follows procedure in Armstrong JQSRT 7, 85 (1967)
            also the same as the intrinsic in the previous version of IDL
    MODIFICATION HISTORY:
            J. Murthy, Mar 1990 (adapted from the FORTRAN program of Armstrong)
                      Sep 1990 (better overflow checking)
    """
    x = np.asarray(u).astype(np.float64)
    y = np.asarray(a).astype(np.float64)

    w = np.array([0.462243670,   0.286675505,   0.109017206,
                  0.0248105209,  0.00324377334, 0.000228338636,
                  7.80255648e-6, 1.08606937e-7, 4.39934099e-10,
                  2.22939365e-13])

    t = np.array([0.245340708, 0.737473729, 1.23407622, 1.73853771,
                  2.25497400,  2.78880606,  3.34785457, 3.94476404,
                  4.60368245,  5.38748089])

    # Hummer's Chebyshev Coefficients
    c = ( 0.1999999999972224, -0.1840000000029998,   0.1558399999965025,
         -0.1216640000043988,  0.0877081599940391,  -0.0585141248086907,
          0.0362157301623914, -0.0208497654398036,   0.0111960116346270,
         -0.56231896167109e-2, 0.26487634172265e-2, -0.11732670757704e-2,
          0.4899519978088e-3, -0.1933630801528e-3,   0.722877446788e-4,
         -0.256555124979e-4,   0.86620736841e-5,    -0.27876379719e-5,
          0.8566873627e-6,    -0.2518433784e-6,      0.709360221e-7,
         -0.191732257e-7,      0.49801256e-8,       -0.12447734e-8,
          0.2997777e-9,       -0.696450e-10,         0.156262e-10,
         -0.33897e-11,         0.7116e-12,          -0.1447e-12,
          0.285e-13,          -0.55e-14,             0.10e-14,
         -0.2e-15)

    y2 = y * y

    # limits are y<1.,  x<4 or y<1.8(x+1),  x>4 (no checking performed)
    u1 = np.exp(-x * x + y2) * np.cos(2. * x * y)

    # Clenshaw's Algorithm
    bno1 = np.zeros(x.shape)
    bno2 = np.zeros(x.shape)
    x1 = np.clip((x / 5.), -np.inf, 1.)
    coef = 4. * x1 * x1 - 2.
    for i in range(33, -1, -1):
        bn = coef * bno1 - bno2 + c[i]
        bno2 = np.copy(bno1)
        bno1 = np.copy(bn)

    f = x1 * (bn - bno2)
    dno1 = 1. - 2. * x * f
    dno2 = f

    q = np.abs(x) > 5
    if q.any():
        x14 = np.power(np.clip(x[q], -np.inf, 500.),  14)
        x12 = np.power(np.clip(x[q], -np.inf, 1000.), 12)
        x10 = np.power(np.clip(x[q], -np.inf, 5000.), 10)
        x8  = np.power(np.clip(x[q], -np.inf, 50000.), 8)
        x6  = np.power(np.clip(x[q], -np.inf, 1.e6),   6)
        x4  = np.power(np.clip(x[q], -np.inf, 1.e9),   4)
        x2  = np.power(np.clip(x[q], -np.inf, 1.e18),  2)
        dno1[q] = -(0.5 / x2 + 0.75 / x4 + 1.875 / x6 +
                    6.5625 / x8 + 29.53125 / x10 +
                    162.4218 / x12 + 1055.7421 / x14)
        dno2[q] = (1. - dno1[q]) / (2. * x[q])

    funct = y * dno1
    if (y > 1.e-8).any():
        q = 1.0
        yn = y
        for i in range(2, 51):
            dn = (x * dno1 + dno2) * (-2. / i)
            dno2 = dno1
            dno1 = dn
            if (i % 2) == 1:
                q = -q
                yn = yn * y2
                g = dn.astype(np.float64) * yn
                funct = funct + q * g
                if np.max(np.abs(g / funct)) <= 1.e-8: break

    k1 = u1 - 1.12837917 * funct
    k1 = k1.astype(np.float64).clip(0)
    return k1

