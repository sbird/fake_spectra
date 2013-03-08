# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

def plot_prochaska_2008_data():
  data = np.loadtxt("vel_width_mtlkin_704.dat")
  vel_data = data[:,2]
  
  nlos = np.size(vel_data)
  
  #These are the same bins Tescari uses
  v_table = np.array([0,40,80,120,160,200,240,320,390,590,800])
  center = np.array([(v_table[i]+v_table[i+1])/2. for i in range(0,np.size(v_table)-1)])
  
  nn = np.histogram(vel_data,v_table)[0]
  #Bin width
  width = np.array([v_table[i+1]-v_table[i] for i in range(0,np.size(v_table)-1)])
  
  #This is the avg. fraction of DLAs per unit absorption distance.
  #It is needed to match the normalisation of Pontzen 2008.
  #DLAdX = 0.065
  
  norm = width* (1.*np.size(vel_data))
  
  vels = nn / norm
  #Use poisson errors
  verr = (np.sqrt(nn))/norm
  #These bins will have error bars that go to 0
  ind = np.where(nn == 1)
  verr[ind] *= (0.98)
  
  plt.errorbar(center,vels,xerr=[center-v_table[:-1],v_table[1:]-center],yerr=verr,fmt='.')

  plt.semilogy(center, vels,'o')
  return (center, vels,verr)
