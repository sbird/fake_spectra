"""Plots testing how the rate network behaves under various approximations."""
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
from rate_network import RateNetwork

def make_cont_plot(zfunc, rates, levels=None, title=""):
    """Helper function to make a contour plot."""
    dens = np.logspace(-6, -2, 100)

    ienergy = np.linspace(100, 700, 100)

    gridd, gridt = np.meshgrid(dens, ienergy)
    temp = rates.get_temp(np.median(gridd)*np.ones_like(gridt), gridt)
    cs = plt.contour(gridd, temp/1e4, zfunc(gridd, gridt), levels=levels)

    plt.clabel(cs,inline=1, fmt="%g")
    plt.xscale('log')

    plt.xlabel(r"$\rho$ (cm$^{-3}$)")
    plt.ylabel(r"T ($10^4$ K)")
    plt.title(title)

def collisplot(zz=3):
    """Plot the effect of collisional ionization on the neutral fraction."""
    rates = RateNetwork(zz,selfshield=False, recomb="C92")
    ratesphoto = RateNetwork(zz,selfshield=False, recomb="C92", collisional=False)

    def ratio(gridd, gridt):
        """Ratio of neutral frac with and without collisional ionizations."""
        nh0 = rates.get_neutral_fraction(gridd, gridt)
        nh0photo = ratesphoto.get_neutral_fraction(gridd, gridt)
        return nh0/nh0photo

    make_cont_plot(ratio, rates, levels=np.array([0.8, 0.9, 0.95, 0.99]), title ="Collisional ionization effect")

def simplenh0frac(gridd, gridt, rates, helium=0.24):
    """Simple formula for nH0/nH neglecting collisions and assuming full ionization of helium."""
    yy = helium/(4-4*helium)
    nh = (1-helium) * gridd
    temp = rates.get_temp(np.median(gridd)*np.ones_like(gridt), gridt)
    phfac = rates.recomb.alphaHp(temp) / rates.photo.gH0(rates.redshift) * (1 + 2 * yy)
    return nh * phfac

def simpleplot(zz=3):
    """Plot the rate network neutral fraction against the simple rescaling formula."""
    rates = RateNetwork(zz,selfshield=False, recomb="C92")

    def ratio(gridd, gridt):
        """Ratio of full neutral frac to simple formula."""
        nh0 = rates.get_neutral_fraction(gridd, gridt)
        nh0simple = simplenh0frac(gridd, gridt, rates)
        return nh0/nh0simple

    make_cont_plot(ratio, rates, levels=np.array([0.8, 0.9, 0.95, 0.99, 1]), title =r"$n_{H0}$ (full) / $n_{H0}$ (simple formula)")

def neplot(zz=3):
    """Plot the rate network electron density against a simple estimate with fully ionized helium."""
    rates = RateNetwork(zz,selfshield=False, recomb="C92")

    def ratio(gridd, gridt, helium=0.24):
        """Ratio of full neutral frac to simple formula."""
        ne = rates.get_ne_by_nh(gridd, gridt)
        yy = helium/(4-4*helium)
        nes = 1 + 2 * yy
        return ne/nes

    make_cont_plot(ratio, rates, title ="ne Simple formula")

def hfracplot(zz=3):
    """Plot the rate network electron density against a simple estimate with fully ionized helium."""
    rates = RateNetwork(zz,selfshield=False, recomb="C92")
    make_cont_plot(rates.get_neutral_fraction, rates, title =r"$n_{H0}/n_{H+}$", levels=np.array([1e-5,1e-4,1e-3,1e-2]))

def hefracplot(zz=3):
    """Plot the rate network electron density against a simple estimate with fully ionized helium."""
    rates = RateNetwork(zz,selfshield=False, recomb="C92")

    def ratio(gridd, gridt, helium=0.24):
        """Ratio of full neutral frac to simple formula."""
        ne = rates.get_equilib_ne(gridd, gridt, helium=helium)
        nh = (1-helium) * gridd
        temp = rates.get_temp(np.median(gridd)*np.ones_like(gridt), gridt)
        nhep = rates._nHep(nh, temp, ne)
        nhepp = rates._nHepp(nh, temp, ne)
        return nhep/nhepp

    make_cont_plot(ratio, rates, title =r"$n_{He+}/n_{He++}$", levels=np.array([0.1, 0.25, 0.5, 1]))

if __name__ == "__main__":
    simpleplot(3)
    plt.savefig("plots/neutral_fraction_simple.pdf")
