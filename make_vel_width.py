#!/usr/bin env python
# -*- coding: utf-8 -*-
"""Make some plots of the velocity widths from the cosmo runs"""

import matplotlib
matplotlib.use('PDF')

import matplotlib.pyplot as plt

import plot_spectra as ps
import vel_data
import leastsq as ls
import os.path as path
import os
import numpy as np
import myname
import math
from save_figure import save_figure

outdir = path.join(myname.base, "plots/")
print "Plots at: ",outdir
zrange = {1:(7,3.5), 3:(3.5,2.5), 5:(2.5,0)}
#Colors and linestyles for the simulations
colors = {0:"pink", 1:"purple", 2:"cyan", 3:"green", 4:"gold", 5:"orange", 7:"blue", 6:"grey", 8:"pink", 9:"red", 'A':"grey"}
colors2 = {0:"darkred", 1:"indigo", 2:"cyan", 3:"darkgreen", 4:"gold", 5:"orange", 7:"darkblue", 6:"grey",8:"cyan", 9:"darkred",'A':"grey"}
lss = {0:"--",1:":", 2:":",3:"-.", 4:"--", 5:"-",6:"--",7:"-", 8:"-",9:"--",'A':"--"}
labels = {0:"ILLUS",1:"HVEL", 2:"HVNOAGN",3:"NOSN", 4:"WMNOAGN", 5:"MVEL",6:"METAL",7:"DEF", 8:"RICH",9:"FAST", 'A':"MOM"}

hspec_cache = {}

def get_hspec(sim, snap, snr=0., box = 25):
    """Get a spectra object, possibly from the cache"""
    halo = myname.get_name(sim, True, box=box)
    #Load from a save file only
    try:
        hspec = hspec_cache[(halo, snap)]
    except KeyError:
        hspec = ps.PlottingSpectra(snap, halo, label=labels[sim], snr=snr)
        hspec_cache[(halo, snap)] = hspec
    return hspec

def plot_vel_width_sim(sim, snap, color="red", HI_cut = None):
    """Load a simulation and plot its velocity width"""
    hspec = get_hspec(sim, snap)
    hspec.plot_vel_width("Si", 2, color=color, HI_cut = HI_cut)

def plot_sep_frac(sim, snap):
    """Plot fraction of lines from separated halos"""
    hspec = get_hspec(sim, snap)
    hspec.plot_sep_frac(color=colors[sim], ls=lss[sim])

def plot_spectrum(sim, snap, num, low=0, high=-1, offset=0,subdir=""):
    """Plot a spectrum"""
    hspec = get_hspec(sim, snap, snr=20., box=10)
    tau = hspec.get_observer_tau("Si", 2, num)
    sdir = path.join(outdir,"spectra/"+subdir)
    if not path.exists(sdir):
        os.mkdir(sdir)
    tau_l = np.roll(tau, offset)
    assert np.max(tau) > 0.1
    hspec.plot_spectrum(tau_l[low:high], flux=False)
    save_figure(path.join(sdir,str(num)+"_cosmo"+str(sim)+"_Si_tau"))
    plt.clf()
#     hspec.plot_spectrum(tau_l[low:high])
#     save_figure(path.join(sdir,str(num)+"_cosmo"+str(sim)+"_Si_spectrum"))
#     plt.clf()
    tau = hspec.get_tau("Si", 2,1260, num)
    tau_l = np.roll(tau, offset)
    hspec.plot_spectrum(tau_l[low:high])
    save_figure(path.join(sdir,str(num)+"_cosmo"+str(sim)+"_Si_1260_spectrum"))
    plt.clf()

def plot_colden(sim, snap, num, subdir="", xlim=100):
    """Plot column density"""
    hspec = get_hspec(sim, snap, snr=20., box=10)
    col_den = hspec.get_col_density("Si", 2)
    mcol = np.max(col_den[num])
    ind_m = np.where(col_den[num] == mcol)[0][0]
    col_den = np.roll(col_den[num], np.size(col_den[num])/2 - ind_m)
    hspec.plot_col_density(col_den)
    plt.ylabel(r"N$_\mathrm{SiII}$ cm$^{-2}$")
    plt.xlim(-1*xlim, xlim)
    plt.ylim(ymin=1e9)
    sdir = path.join(outdir,"spectra/"+subdir)
    if not path.exists(sdir):
        os.mkdir(sdir)
    save_figure(path.join(sdir,str(num)+"_cosmo"+str(sim)+"_Si_colden"))
    plt.clf()
    col_den = hspec.get_col_density("H", 1)
    col_den = np.roll(col_den[num], np.size(col_den[num])/2 - ind_m)
    hspec.plot_col_density(col_den)
    plt.xlim(-1*xlim, xlim)
    plt.ylabel(r"N$_\mathrm{HI}$ cm$^{-2}$")
    plt.ylim(ymin=1e15)
    save_figure(path.join(sdir,str(num)+"_cosmo"+str(sim)+"_H_colden"))
    plt.clf()

def plot_spectrum_max(sim, snap, velbin, velwidth, num):
    """Plot spectrum with max vel width"""
    hspec = get_hspec(sim, snap, snr=20., box=10)
    vels = hspec.vel_width("Si",2)
    ind = hspec.get_filt("Si", 2)
    subdir = path.join("cosmo"+str(sim)+"-10",str(velbin))
    band = np.intersect1d(ind[0], np.where(np.logical_and(vels > velbin-velwidth, vels < velbin+velwidth))[0])
    np.random.seed(2323)
    index = np.random.randint(0, np.size(band), num)
#     (low, high, offset) = hspec.find_absorber_width("Si",2, minwidth=1.2*velbin/2.)
    (low, high, offset) = hspec.find_absorber_width("Si",2, minwidth=1.1*(velbin+velwidth))
    for nn in band[index]:
        plot_spectrum(sim, snap, nn, low[nn], high[nn], offset[nn], subdir=subdir)
        plot_colden(sim, snap, nn, subdir, xlim=np.max((vels[nn]/2.,100)))

def plot_metallicity(sims, snap):
    """Plot metallicity, vel width, their correlation and the extra statistics"""
    out = "cosmo_metallicity_z"+str(snap)
    for sim in sims:
        hspec = get_hspec(sim, snap)
        hspec.plot_metallicity(color=colors[sim], ls=lss[sim])
    vel_data.plot_alpha_metal_data(zrange[snap])
    plt.legend(loc=2,ncol=3)
    plt.ylim(0,2)
    save_figure(path.join(outdir,out))
    plt.clf()

def plot_met_corr(sims,snap):
    """Plot metallicity vel width correlations"""
    for sim in sims:
        out = "cosmo"+str(sim)+"_correlation_z"+str(snap)
        hspec = get_hspec(sim, snap)
        hspec.plot_Z_vs_vel_width(color=colors[sim], color2=colors2[sim])
        vel_data.plot_prochaska_2008_correlation(zrange[snap])
        plt.xlim(10, 500)
        save_figure(path.join(outdir,out))
        plt.clf()

def plot_vel_width_sims(sims, snap, log=False):
    """Plot velocity widths for a series of simulations"""
    vel_data.plot_prochaska_2008_data()
    for sss in sims:
        #Make abs. plot
        hspec = get_hspec(sss, snap)
        hspec.plot_vel_width("Si", 2, color=colors[sss], ls=lss[sss])
    outstr = "cosmo_vel_width_z"+str(snap)
    if log:
        ax = plt.gca()
        ax.set_yscale('log')
        plt.ylim(1e-2,10)
        outstr+="_log"
    else:
        plt.ylim(1e-2,2)
    plt.xlabel(r"$v_\mathrm{90}$ (km s$^{-1}$)")
    plt.xlim(10,1000)
    plt.legend(loc=2,ncol=3)
    save_figure(path.join(outdir,outstr))
    plt.clf()

def plot_eq_width(sims, snap):
    """Plot velocity widths for a series of simulations"""
    for sss in sims:
        #Make abs. plot
        hspec = get_hspec(sss, snap)
        hspec.plot_eq_width("Si", 2, 1526, color=colors[sss], ls=lss[sss])
    outstr = "cosmo_eq_width_z"+str(snap)
    vel_data.plot_si1526_eqw(zrange[snap], nv_table=7)
    plt.ylim(0,3)
    plt.legend(loc=2,ncol=3)
    save_figure(path.join(outdir,outstr))
    plt.clf()

def plot_rel_vel_width(sims, snap):
    """Plot velocity widths relative to simulation 7"""
    hspec = get_hspec(7, snap)
    (vbin, vels7) = hspec.vel_width_hist("Si", 2)
    #Make rel plot
    for sss in sims:
        hspec = get_hspec(sss, snap)
        (vbin, vel) = hspec.vel_width_hist("Si", 2)
        mm = np.min([np.size(vel), np.size(vels7)])
        plt.semilogx(vbin[:mm], vel[:mm]/vels7[:mm], color=colors[sss],ls=lss[sss])
    plt.xlim(10, 1000)
    save_figure(path.join(outdir,"cosmo_rel_vel_z"+str(snap)))
    plt.clf()

def plot_mean_median(sims, snap):
    """Plot mean-median statistic for all sims on one plot"""
    #Plot extra statistics
    for sss in sims:
        hspec = get_hspec(sss, snap)
        hspec.plot_f_meanmedian("Si", 2, color=colors[sss], ls=lss[sss])
    vel_data.plot_extra_stat_hist(False)
    plt.ylim(0,3)
    plt.legend(loc=2,ncol=3)
    save_figure(path.join(outdir,"cosmo_mean_median_z"+str(snap)))
    plt.clf()

import numexpr as ne

class RotationFiltered(ps.PlottingSpectra):
    """Class to plot the velocity widths of only rotationally supported gas"""
    def _filter_particles(self, elem_den, pos, velocity, den):
        """Filtered list of particles that are rotationally supported by a halo."""
        #Filter particles that are non-dense, as they will not be in halos
        ind2 = np.where(np.logical_and(den > 3e-4, elem_den > 0))
        halo_cofm = self.sub_cofm
        sub_cofm = self.sub_sub_cofm
        ind3 = []
        frachigh = 1.5
        fraclow = 0.7
        non_rot = 0
        for ii in ind2[0]:
            #Is this within the virial radius of any halo?
            ppos = pos[ii,:]
            dd = ne.evaluate("sum((halo_cofm - ppos)**2,axis=1)")
            ind = np.where(dd < self.sub_radii**2)
            #Check subhalos
            if np.size(ind) < 1:
                dd = ne.evaluate("sum((sub_cofm - ppos)**2,axis=1)")
                ind = np.where(dd < self.sub_sub_radii**2)
                if np.size(ind) < 1:
                    continue
                ind = ind[0][0]
                hvel = self.sub_sub_vel[ind,:]
                hcofm = self.sub_sub_cofm[ind,:]
                hrad = self.sub_sub_radii[ind]
                vvir = self.virial_vel([ind,], subhalo=True)[0]
            else:
                ind = ind[0][0]
                hvel = self.sub_vel[ind,:]
                hcofm = self.sub_cofm[ind,:]
                hrad = self.sub_radii[ind]
                vvir = self.virial_vel([ind,])[0]
            #It is! What is the perpendicular velocity wrt this halo?
            lvel =  velocity[ii, :] - hvel
            #Radial vector from halo
            lpos = ppos - hcofm
            ldist = np.sqrt(np.sum(lpos**2))
            #Find parallel by dotting with unit vector
            vpar = np.dot(lvel, lpos/ldist)
            vperp = np.sqrt(np.sum(lvel**2) - vpar)
            #Rotational velocity assuming NFW concentration 10 (like MW).
            vhalo = vvir * np.sqrt(5*ldist) / (1+ 10 * ldist / hrad)
            #Are we rotation supported?
            #Also, the angular vector should dominate over the radial
            if np.abs(vpar / vperp) > 0.3 or vperp / vhalo >  frachigh or vperp / vhalo < fraclow:
                non_rot += 1
                continue
            #If we are, add to the list
            ind3+= [ii,]
        print "Filtered ",np.size(ind2[0])," particles to ",np.size(ind3)
        print "Non-rotating ",non_rot
        return ind3

    def get_filt(self, elem, ion, thresh = 2e11):
        """
        Get an index list to exclude spectra where the ion is not observable.
        For DLAs this should not be a huge fraction of the total spectra,
        as few metal-poor DLAs are observed.

        thresh - observable column density threshold
        """
        #Remember this is not in log...
        #This is in practice for SiII ~200/5000 ~ 4% of spectra, which is a little large.
        #Oddly it is less for metallicity
        met = np.max(self.get_col_density(elem, ion), axis=1)
        ind = np.where(np.logical_and(met > thresh, np.max(self.get_observer_tau(elem, ion), axis=1) > 0.1))
        print "Sightlines with rotating absorption: ",np.size(ind)
        return ind

def plot_v_struct(sims, snap):
    """Plot mean-median statistic for all sims on one plot"""
    #Plot extra statistics
    for sss in sims:
        halo = myname.get_name(sss, True)
        #Load from a save file only
        hspec = RotationFiltered(snap, halo, label=labels[sss], spec_res = 0.1)
        hspec.get_observer_tau("Si",2,force_recompute=True)
        hspec.plot_vel_width("Si", 2, color=colors[sss], ls=lss[sss])
    vel_data.plot_prochaska_2008_data()
    plt.legend(loc=2,ncol=3)
    plt.xlabel(r"$v_\mathrm{90}$ km s$^{-1}$")
    plt.xlim(2, 1000)
    plt.ylim(0,2)
    save_figure(path.join(outdir, "cosmo_rot_z"+str(snap)))
    plt.clf()

def plot_f_peak(sims, snap):
    """Plot peak statistic for all sims on one plot"""
    for sss in sims:
        hspec = get_hspec(sss, snap)
        hspec.plot_f_peak("Si", 2, color=colors[sss], ls=lss[sss])
    plt.legend(loc=2,ncol=3)
    vel_data.plot_extra_stat_hist(True)
    plt.ylim(0,3)
    save_figure(path.join(outdir,"cosmo_peak_z"+str(snap)))
    plt.clf()

def plot_vel_widths_cloudy():
    """Plot some velocity width data for different cloudy models"""
    #Load sims
    hspec0 = ps.PlottingSpectra(3, myname.get_name(0, True))
    hspec1 = ps.PlottingSpectra(3, myname.get_name(0,True), savefile="rand_spectra_DLA_fancy_atten.hdf5")
    #Make abs. plot
    hspec0.plot_vel_width("Si", 2, color="blue", ls="--")
    hspec1.plot_vel_width("Si", 2, color="red", ls="-")
    vel_data.plot_prochaska_2008_data()
    save_figure(path.join(outdir,"cosmo_feedback_cloudy_z3"))
    plt.clf()
    #Make rel plot
    (vbin, vels0) = hspec0.vel_width_hist("Si", 2)
    (vbin, vels2) = hspec1.vel_width_hist("Si", 2)
    mm = np.min((np.size(vels2),np.size(vels0)))
    plt.semilogx(vbin[:mm], vels0[:mm]/vels2[:mm], color="blue",ls="-")
    plt.xlim(1, 1000)
    save_figure(path.join(outdir,"cosmo_rel_vel_cloudy_z3"))
    plt.clf()

def plot_vel_redshift_evo(sim):
    """Plot the evolution with redshift of a simulation"""
    halo = myname.get_name(sim, True)
    vels = {}
    for snap in (1,3,5):
        hspec0 = ps.PlottingSpectra(snap, halo)
        (vbin, vels[snap]) = hspec0.vel_width_hist("Si", 2)
    mm = np.min([np.size(vel) for vel in vels.values()])
    #Normalised by z=3
    plt.semilogx(vbin[:mm], vels[5][:mm]/vels[3][:mm], color="black",ls="--")
    plt.semilogx(vbin[:mm], vels[1][:mm]/vels[3][:mm], color="grey",ls="-")
    plt.xlim(10, 1000)
    plt.ylim(0.5,1.5)
    save_figure(path.join(outdir,"cosmo"+str(sim)+"_zz_evol"))
    plt.clf()

def disc_distrib(phi, iota, scale=0.1):
    """
    Get the velocity width from an iota and phi.
    Arguments: phi - azimuthal angle
              iota - inclination angle
              scale - ratio of disc scale length to scale height, h / R_d.
    """
    return 2 * np.sin(phi) / (1 + np.tan(iota) / np.cos(phi) / scale)

def make_disc_model(nsample = 50000, scale=0.1,color="red", label="Disc"):
    """Plot the distribution of velocity width by virial velocity that would result from a rotating disc"""
    iota = math.pi/2.*np.random.random_sample(nsample)
    phi = 2.*math.pi*np.random.random_sample(nsample)
    v90 = disc_distrib(phi, iota,scale)
    v_table = 10**np.arange(-3, np.log10(np.max(v90)), 0.1)
    vels = np.histogram(np.log10(v90), np.log10(v_table),density=True)[0]
    vbin = np.array([(v_table[i]+v_table[i+1])/2. for i in range(0,np.size(v_table)-1)])
    plt.semilogx(vbin, vels, color=color, lw=3, ls="--",label=label)

def read_H_model():
    """Read and plot the data from Haehnelt et al 1998"""
    data = np.loadtxt("damp12_f.dat")
    v90 = data[:,2]
    vvir = data[:,10]
    v_table = 10**np.arange(-2, np.log10(np.max(v90/vvir)), 0.1)
    vels = np.histogram(np.log10(v90/vvir), np.log10(v_table),density=True)[0]
    vbin = np.array([(v_table[i]+v_table[i+1])/2. for i in range(0,np.size(v_table)-1)])
    plt.semilogx(vbin, vels, color="purple", lw=3, ls=":",label="Haehnelt 98")

def plot_vvir_models():
    """Plot histogram of velocity width by virial velocity"""
    #Load from a save file only
    hspec = get_hspec(7,3)
    hspec.plot_virial_vel_vs_vel_width("Si", 2, color=colors[7], ls=lss[7], label=labels[7])
    hspec = get_hspec(3,3)
    hspec.plot_virial_vel_vs_vel_width("Si", 2, color=colors[3], ls=lss[3], label=labels[3])
    make_disc_model(scale=0.25,label="Disc")
    read_H_model()
    plt.legend(loc=2)
    plt.xlim(0.01, 10)
    plt.xlabel(r"$v_\mathrm{90} / v_\mathrm{vir}$")
    save_figure(path.join(outdir, "vvir90_model"))

def do_statistics(sim, snap):
    """Compute statistics"""
    #Get Observational data
    (_, met, vel) = vel_data.load_data(zrange[snap])
    vel = np.log10(vel)
    #Get Simulated data
    halo = myname.get_name(sim, True)
    hspec = ps.PlottingSpectra(snap, halo)
    svel = hspec.vel_width("Si", 2)
    smet = hspec.get_metallicity()
    #Ignore objects too faint to be seen
    ind2 = np.where(smet > 1e-4)
    smet = np.log10(smet[ind2])
    svel = np.log10(svel[ind2])
    #Fit to both datasets
    (obs_intercept, obs_slope, obs_var) = ls.leastsq(vel,met)
    (s_intercept, s_slope, s_var) = ls.leastsq(svel,smet)
    print "obs fit: ",obs_intercept, obs_slope, np.sqrt(obs_var)
    print "sim fit: ",s_intercept, s_slope, np.sqrt(s_var)
    #Find correlations
    print "obs pearson r: ",ls.pearson(vel, met,obs_intercept, obs_slope)
    print "sim pearson r: ",ls.pearson(svel, smet,s_intercept, s_slope)
    print "obs kstest: ",ls.kstest(vel, met,obs_intercept, obs_slope)
    print "sim kstest: ",ls.kstest(svel, smet,s_intercept, s_slope)
    #Now test whether they come from the same population
    kss = hspec.kstest(10**met, 10**vel)
    print "KS test between simulated and observed samples: ",kss
    #Do 200 trials and see how many times the KS test is worse
    ntrials = 50
    count = 0
    for _ in xrange(ntrials):
        rand = np.random.randint(0,np.size(svel), np.size(vel))
        if kss <= hspec.kstest(10**smet[rand], 10**svel[rand]):
            count+=1
    print "Prob KS test between simulated samples was larger: ",count*1./ntrials

if __name__ == "__main__":
#     plot_vel_widths_cloudy()

    for zz in (1,3,5):
        do_statistics(7,zz)
    for ss in (1,3,9):
        do_statistics(ss,3)

    plot_spectrum_max(5,3, 60, 20, 15)
    plot_spectrum_max(5,3, 100, 20, 15)
    plot_spectrum_max(5,3, 200, 35, 15)
    plot_spectrum_max(5,3, 400, 50, 15)
    simlist = (1,3,7,9) #range(8)

#     for ss in simlist:
#         plot_spectrum_max(ss,3, 60, 20, 15)
#         plot_spectrum_max(ss,3, 200, 35, 15)
#         plot_spectrum_max(ss,3, 400, 50, 15)

#     plot_spectrum_max(7,3, 100, 20, 15)
#     plot_spectrum_max(7,3, 140, 20, 15)
    plot_vel_width_sims(simlist, 3, log=True)
    for zz in (1,3,5):
        plot_met_corr(simlist,zz)
        hspec_cache = {}

    for zz in (1, 3, 5):
#         plot_v_struct(simlist, zz)
        plot_eq_width(simlist, zz)
        plot_metallicity(simlist, zz)
        plot_vel_width_sims(simlist, zz)
        plot_mean_median(simlist, zz)
        plot_f_peak(simlist, zz)

    for ss in simlist:
        plot_sep_frac(ss,3)
    plt.legend(loc=2,ncol=3)
    save_figure(path.join(outdir,"cosmo_sep_frac_z3"))
    plt.clf()

    for ss in simlist:
        plot_vel_redshift_evo(ss)
