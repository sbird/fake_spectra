#!/usr/bin/env python2
"""Short script to read cloudy table files and output numpy arrays."""

import numpy as np
import re
import os.path as path
from scipy.interpolate import interp1d

#Max number of ion species to count
nions = 17

def handle_single_line(splitline):
    """Short function to handle a single cloudy line, split by whitespace, and deal with lack of spaces"""
    #Sometimes Cloudy fails to put a space between output columns.
    #Test for this by trying to create a float from the element
    #and seeing if it fails.
    ii = 0
    while ii < np.size(splitline):
        s = splitline[ii]
        try:
            s=float(s)
        except ValueError:
            #Cloudy elements are <= 0, so we can split on -.
            tmp = s.split("-")
            splitline[ii] = -1*float(tmp[-1])
            #Insert elements backwards: First element is ignored because it will be ""
            for q in range(-2,-np.size(tmp),-1):
                splitline.insert(ii, -1*float(tmp[q]))
                ii+=1
        ii+=1
    return splitline

def convert_single_file(ion_table):
    """Convert a single file to a numpy (text) format.
    Discard all metals except the 8 we follow in Arepo, which are:
    He, C, N, O, Ne, Mg, Si, Fe
    Output an 8 x 22 array, the number of ions of iron cloudy gives.
    Array elements greater than the physical ionisation of the metal will be set to -30
    which is what cloudy uses for log(0).
    """
    #Species names to search the cloudy table for
    species = ("Helium", "Carbon", "Nitrogen", "Oxygen", "Neon", "Magnesium", "Silicon", "Iron")
    #Next element
    elmt = 0
    num_table = -30*np.ones((np.size(species), nions))

    f=open(ion_table)
    line = f.readline()
    while line != "" and elmt < np.size(species):
        #This line contains a species we want
        if re.search(species[elmt],line):
            #Split with whitespace
            splitline=line.split()
            #Remove name column
            if splitline[0] == species[elmt]:
                splitline = splitline[1:]
            else:
                splitline[0] = re.sub(species[elmt], "", splitline[0])
            splitline = handle_single_line(splitline)
            #Set table
            num_table[elmt, 0:np.size(splitline)] = np.array(splitline)
            #Iron will continue onto a second line.
            #Ignore ion species > 17. This only affects iron at high densities
            #if species[elmt] == "Iron":
            #    line2 = f.readline()
            #    splitline2 = line2.split()
            #    splitline2 = handle_single_line(splitline2)
            #    start = np.size(splitline)
            #    num_table[elmt, start:(start+np.size(splitline2))] = np.array(splitline2)
            elmt+=1
        line = f.readline()
    return num_table

def read_all_tables(nred, nmet, nrho, directory):
    """
    Read a batch of tables on a grid of redshift (nred), metallicity (nmet) and density (nrho) from a directory.
    Subdirs are in the form: $RED_$MET_$DENS/
    Returns a table with the indices:
    SPECIES, REDSHIFT, METALLICITY, DENSITY, ION
    """
    #Last two indices are num. species and n. ions
    tables = np.empty([8, nred, nmet, nrho,nions])
    for zz in range(1,nred+1):
        for ZZ in range(1, nmet+1):
            for rr in range(1, nrho+1):
                strnums = str(zz)+"_"+str(ZZ)+"_"+str(rr)
                ionfile = directory+"/"+strnums+"/ionization_"+strnums+".dat"
                tables[:,zz-1,ZZ-1,rr-1,:] = convert_single_file(ionfile)
    return tables

class CloudyTable:
    """Class to interpolate tables from cloudy for a desired density metallicity and redshift"""
    def __init__(self,directory):
        """Read a cloudy table from somewhere"""
        self.savefile = path.join(directory,"cloudy_table.npz")
        self.reds = np.array([1,2,3,4,5,6])
        self.mets = np.array([-4,-3,-2,-1,0,1])
        self.dens = np.array([-6,-5,-4,-3,-2,-1,0,1,2])
        self.species = ("He", "C", "N", "O", "Ne", "Mg", "Si", "Fe")
        #Solar abundances from Hazy table 7.1 as Z = n/n(H)
        self.solar = np.array([0.1, 2.45e-4, 8.51e-5, 4.9e-4, 1.e-4, 3.47e-5, 3.47e-5, 2.82e-5])
        try:
            datafile = np.load(self.savefile)
            self.table = datafile["table"]
        except (IOError, KeyError):
            self.table = read_all_tables(np.size(self.reds), np.size(self.mets), np.size(self.dens), directory)
            self.save_file()
        self.directory = directory

    def ion(self,species, ion, red, met, rho):
        """Interpolate a table onto given redshift, metallicity and density for species
        Returns a log( ionisation fraction ).
        rho is the density of hydrogen in atoms / cm^3. Internally the log is taken
        met is the mass-weighted metallicity of a species
        specified as Z = n/ n(H), and internally converted into
        X = log_10(Z / Z_solar) for cloudy
        red is the redshift
        species is the element name of the metal species we want
        ion is the ionisation number
        Returns the fraction of this ion in the species
        """
        cspe = self.species.index(species)
        ind = np.where(np.logical_and(rho != 0, met != 0))
        crho = np.log10(rho[ind])
        cmet = np.log10(met[ind]/self.solar[cspe])
        red_ints = interp1d(self.reds, self.table[cspe,:,:,:,ion],axis = 0)
        met_ints = interp1d(self.mets, red_ints(red),axis = 0)
        #For super small metallicities
        #use the lowest ion fraction we have
        if cmet < np.min(self.mets):
            cmet = np.min(self.mets)
        dens_ints = interp1d(self.dens, met_ints(cmet),axis = 0)
        #For super small densities it won't make a difference,
        #so use the lowest ion fraction we have
        if crho < np.min(self.dens):
            crho = np.min(self.dens)
        ions = np.zeros(np.shape(rho))
        ions[ind] = dens_ints(crho)
        return 10**ions

    def save_file(self):
        """
        Saves tables to a file, because they are slow to generate.
        File is hard-coded to be directory/cloudy_table.npz.
        """
        np.savez(self.savefile,table=self.table)



