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
import numpy as np
import myname
from save_figure import save_figure

outdir = path.join(myname.base, "plots/")
print "Plots at: ",outdir
zrange = {1:(7,3.5), 3:(3.5,2.5), 5:(2.5,0)}
#Colors and linestyles for the simulations
colors = {0:"pink", 1:"purple", 2:"cyan", 3:"green", 4:"gold", 5:"orange", 7:"blue", 6:"grey", 8:"pink", 9:"red", 'A':"grey"}
colors2 = {0:"darkred", 1:"indigo", 2:"cyan", 3:"darkgreen", 4:"gold", 5:"orange", 7:"darkblue", 6:"grey",8:"cyan", 9:"darkred",'A':"grey"}
lss = {0:"--",1:":", 2:":",3:"-.", 4:"--", 5:"-",6:"--",7:"-", 8:"-",9:"--",'A':"--"}
labels = {0:"ILLUS",1:"HVEL", 2:"HVNOAGN",3:"NOSN", 4:"WMNOAGN", 5:"MVEL",6:"METAL",7:"2xUV", 8:"RICH",9:"FAST", 'A':"MOM"}

hspec_cache = {}

def get_hspec(sim, snap):
    """Get a spectra object, possibly from the cache"""
    halo = myname.get_name(sim, True)
    #Load from a save file only
    try:
        hspec = hspec_cache[(sim, snap)]
    except KeyError:
        hspec = ps.PlottingSpectra(snap, halo, label=labels[sim])
        hspec_cache[(sim, snap)] = hspec
    return hspec

def plot_vel_width_sim(sim, snap, color="red", HI_cut = None):
    """Load a simulation and plot its velocity width"""
    hspec = get_hspec(sim, snap)
    hspec.plot_vel_width("Si", 2, color=color, HI_cut = HI_cut)

def plot_sep_frac(sim, snap):
    """Plot fraction of lines from separated halos"""
    hspec = get_hspec(sim, snap)
    hspec.plot_sep_frac(color=colors[sim], ls=lss[sim])

def plot_spectrum(sim, snap, num, subdir="", minwidth=0):
    """Plot a spectrum"""
    hspec = get_hspec(sim, snap)
    tau = hspec.get_observer_tau("Si", 2)
    ind = hspec.get_filt("Si", 2)
    (low,high, offset) = hspec.find_absorber_width("Si",2, minwidth=minwidth)
    tau_l = np.roll(tau[ind][num], offset[ind][num])
    hspec.plot_spectrum(tau_l[low[ind][num]:high[ind][num]], flux=False)
    save_figure(path.join(outdir,"spectra/"+subdir+str(num)+"_cosmo"+str(sim)+"_Si_tau"))
    plt.clf()
    hspec.plot_spectrum(tau_l[low[ind][num]:high[ind][num]])
    save_figure(path.join(outdir,"spectra/"+subdir+str(num)+"_cosmo"+str(sim)+"_Si_spectrum"))
    plt.clf()
    tau = hspec.get_tau("Si", 2,1260)
    tau_l = np.roll(tau[ind][num], offset[ind][num])
    hspec.plot_spectrum(tau_l[low[ind][num]:high[ind][num]])
    save_figure(path.join(outdir,"spectra/"+subdir+str(num)+"_cosmo"+str(sim)+"_Si_1260_spectrum"))
    plt.clf()

def plot_colden(sim, snap, num, subdir="", xlim=100):
    """Plot column density"""
    hspec = get_hspec(sim, snap)
    ind = hspec.get_filt("Si", 2)
    col_den = hspec.get_col_density("Si", 2)[ind]
    mcol = np.max(col_den[num])
    ind_m = np.where(col_den[num] == mcol)[0][0]
    col_den = np.roll(col_den[num], np.size(col_den[num])/2 - ind_m)
    hspec.plot_col_density(col_den)
    plt.ylabel(r"N$_\mathrm{SiII}$ cm$^{-2}$")
    plt.xlim(-1*xlim, xlim)
    plt.ylim(ymin=1e9)
    save_figure(path.join(outdir,"spectra/"+subdir+str(num)+"_cosmo"+str(sim)+"_Si_colden"))
    plt.clf()
    col_den = hspec.get_col_density("H", 1)[ind]
    col_den = np.roll(col_den[num], np.size(col_den[num])/2 - ind_m)
    hspec.plot_col_density(col_den)
    plt.xlim(-1*xlim, xlim)
    plt.ylabel(r"N$_\mathrm{HI}$ cm$^{-2}$")
    plt.ylim(ymin=1e15)
    save_figure(path.join(outdir,"spectra/"+subdir+str(num)+"_cosmo"+str(sim)+"_H_colden"))
    plt.clf()

def plot_spectrum_max(sim, snap):
    """Plot spectrum with max vel width"""
    hspec = get_hspec(sim, snap)
    vels = hspec.vel_width("Si",2)
    ind = hspec.get_filt("Si", 2)
    subdir = "max/cosmo"+str(sim)+"/"
    num = np.where(vels[ind] == np.max(vels[ind]))[0][0]
    plot_spectrum(sim, snap, num, subdir, minwidth=500)
    plot_colden(sim, snap, num, subdir, 500)
    num = np.where(vels[ind] > 450)[0]
    for nn in num:
        plot_spectrum(sim, snap, nn, subdir, minwidth=500)
        plot_colden(sim, snap, nn, subdir,500)


def plot_spectrum_density_velocity(sim, snap, num):
    """Plot a spectrum"""
    hspec = get_hspec(sim, snap)
    hspec.plot_spectrum_density_velocity("Si",2, num)
    save_figure(path.join(outdir,"spectra/cosmo"+str(sim)+"_tdv_Si_spectrum"))
    plt.clf()

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

def plot_v_struct(sims, snap):
    """Plot mean-median statistic for all sims on one plot"""
    #Plot extra statistics
    for sss in sims:
        hspec = get_hspec(sss, snap)
        plt.figure(1)
        hspec.plot_velocity_amp("H",1, color=colors[sss], ls=lss[sss])
        plt.figure(2)
        hspec.plot_velocity_theta("H",1, color=colors[sss], ls=lss[sss])
    plt.figure(1)
    plt.legend(loc=2,ncol=3)
    save_figure(path.join(outdir, "cosmo_amp_z"+str(snap)))
    plt.clf()
    plt.figure(2)
    plt.legend(loc=2,ncol=3)
    save_figure(path.join(outdir, "cosmo_theta_z"+str(snap)))
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
    for i in xrange(ntrials):
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

    simlist = (1,3,7,9) #range(8)
    for ss in simlist:
        for nn in (272,350,457,1030,1496,2030,3333):
            plot_spectrum(ss,3, nn, subdir = "cosmo"+str(ss)+"/")
            plot_colden(ss,3,nn)
        plot_spectrum_max(ss,3)

    plot_vel_width_sims(simlist, 3, log=True)
    for zz in (1,3,5):
        plot_met_corr(simlist,zz)
        hspec_cache = {}

    for zz in (1, 3, 5):
        plot_v_struct(simlist, zz)
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
