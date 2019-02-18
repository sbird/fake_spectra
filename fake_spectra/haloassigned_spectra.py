"""This module contains a class and functions which do analysis on spectra associated to
galactic halos. This is fundamentally a somewhat wooly idea, because absorbers and halos are not always closely associated!"""
from __future__ import print_function
import numpy as np
try:
    import numexpr as ne
except ImportError:
    #Non-essential
    pass
import matplotlib.pyplot as plt

from . import subfindhdf
from . import spec_utils
from . import plot_spectra as ps

try:
    xrange(1)
except NameError:
    xrange = range

class HaloAssignedSpectra(ps.PlottingSpectra):
    """Class which extends the Spectra class to include methods that connect each absorber in a sightline to a galactic halo."""
    def __init__(self, *args, **kwargs):
        ps.PlottingSpectra.__init__(self, *args, **kwargs)
        #Try to load a halo catalogue
        self.load_halo()

    def mass_hist(self, dm=0.1):
        """
        Compute a histogram of the host halo mass of each DLA spectrum.

        Parameters:
            dm - bin spacing

        Returns:
            (mbins, pdf) - Mass (binned in log) and corresponding PDF.
        """
        (halos, _) = self.find_nearest_halo()
        f_ind = np.where(halos != -1)
        #nlos = np.shape(vel_width)[0]
        #print('nlos = ',nlos)
        virial = self.virial_vel(halos[f_ind])
        m_table = 10**np.arange(np.log10(np.min(virial)+0.1), np.log10(np.max(virial)), dm)
        mbin = np.array([(m_table[i]+m_table[i+1])/2. for i in range(0,np.size(m_table)-1)])
        pdf = np.histogram(np.log10(virial),np.log10(m_table), density=True)[0]
        print("Field DLAs: ",np.size(halos)-np.size(f_ind))
        return (mbin, pdf)

    def virial_vel(self, halos=None, subhalo=False):
        """Get the virial velocities of the selected halos in km/s"""
        if subhalo:
            if halos is not None:
                mm = self.sub_sub_mass[halos]
                rr = np.array(self.sub_sub_radii[halos])
            else:
                mm = self.sub_sub_mass
                rr = np.array(self.sub_sub_radii)
        else:
            if halos is not None:
                mm = self.sub_mass[halos]
                rr = np.array(self.sub_radii[halos])
            else:
                mm = self.sub_mass
                rr = np.array(self.sub_radii)
        #physical cm from comoving kpc/h
        cminkpch = self.units.UnitLength_in_cm/self.hubble/(1+self.red)
        # Conversion factor from M_sun/kpc/h to g/cm
        conv = self.units.UnitMass_in_g/1e10/cminkpch
        #Units: grav is in cm^3 /g/s^-2
        #Define zero radius, zero mass halos as having zero virial velocity.
        rr[np.where(rr == 0)] = 1
        virial = np.sqrt(self.units.gravcgs*conv*mm/rr)/1e5
        return virial

    def min_halo_mass(self, minpart = 400):
        """Min resolved halo mass in internal Gadget units (1e10 M_sun)"""
        #This is rho_c in units of h^-1 1e10 M_sun (kpc/h)^-3
        rhom = 2.78e+11* self.OmegaM / 1e10 / (1e3**3)
        #Mass of an SPH particle, in units of 1e10 M_sun, x omega_m/ omega_b.
        try:
            target_mass = self.box**3 * rhom / self.npart[0]
        except AttributeError:
            #Back-compat hack
            target_mass = self.box**3 * rhom / 512.**3
        min_mass = target_mass * minpart
        return min_mass

    def load_halo(self):
        """Load a halo catalogue: note this will return some halos with zero radius.
           These are empty FoF groups and should not be a problem."""
        SolarMass_in_g=1.989e33
        try:
            subs=subfindhdf.SubFindHDF5(self.base, self.num)
            #Get particle center of mass, use group catalogue.
            self.sub_cofm=subs.get_grp("GroupPos")
            #halo masses in M_sun/h: use M_200
            self.sub_mass=subs.get_grp("Group_M_Crit200")*self.units.UnitMass_in_g/SolarMass_in_g
            #r200 in kpc/h (comoving).
            self.sub_radii = subs.get_grp("Group_R_Crit200")
            self.sub_vel = subs.get_grp("GroupVel")
            self.sub_sub_radii =  subs.get_sub("SubhaloHalfmassRad")
            self.sub_sub_cofm =  subs.get_sub("SubhaloPos")
            self.sub_sub_mass =  subs.get_sub("SubhaloMass")
            self.sub_sub_index = subs.get_sub("SubhaloGrNr")
            self.sub_sub_vel = subs.get_sub("SubhaloVel")
        except IOError:
            pass

    def assign_to_halo(self, zpos, halo_radii, halo_cofm):
        """
        Assign a list of lists of positions to halos, by finding the unique halo
        within whose virial radius each position is.
        """
        dists = []
        halos = []
        #X axis first
        for ii in xrange(len(zpos)):
            proj_pos = np.array(self.cofm[ii,:])
            ax = self.axis[ii]-1
            dists.append([])
            halos.append([])
            for zzp in zpos[ii]:
                proj_pos[ax] = zzp
                #Is this within the virial radius of any halo?
                try:
                    dd = ne.evaluate("sum((halo_cofm - proj_pos)**2,axis=1)")
                except NameError:
                    dd = np.sum((halo_cofm - proj_pos)**2,axis=1)
                ind = np.where(dd < halo_radii**2)
                #Should not be multiple close halos
                # assert(np.size(ind) < 2)
                #Very rarely, in 2/5000 cases,
                #something hits the edge of more than one halo
                #This is so rare we don't worry about it.
                if np.size(ind) >= 1:
                    halos[ii].append(ind[0][0])
                    dists[ii].append(np.sqrt(dd[ind][0]))
        return (halos, dists)

    def get_contiguous_regions(self, elem="H", ion = 1, thresh = 2e20, relthresh = 1e-3):
        """
        Find the weighted z position of all contiguous DLA-hosting regions in each spectrum.
        Returns a list of lists. Each element in the outer list corresponds to a spectrum.
        Each inner list is the list of weighted z positions of DLA-hosting regions.
        """
        den = self.get_col_density(elem, ion)
        contig = []
        seps = np.zeros(self.NumLos, dtype=np.bool)
        (roll, colden) = spec_utils.get_rolled_spectra(den)
        #deal with periodicity by making sure the deepest point is in the middle
        for ii in xrange(self.NumLos):
            # This is column density, not absorption, so we cannot
            # use the line width to find the peak region.
            lcolden = colden[ii,:]
            # Get first and last indices of separate regions in list
            if np.max(lcolden) > thresh:
                seps = combine_regions(lcolden > thresh)
            else:
                seps = combine_regions(lcolden > relthresh*np.max(lcolden))
            # Find weighted z position for each one
            zposes = []
            for jj in xrange(np.shape(seps)[0]):
                nn = np.arange(self.nbins)[seps[jj,0]:seps[jj,1]]-roll[ii]
                llcolden = lcolden[seps[jj,0]:seps[jj,1]]
                zpos = np.sum(llcolden*nn)
                summ = np.sum(llcolden)
                #Make sure it refers to a valid position
                zpos = (zpos / summ) % self.nbins
                zpos *= 1.*self.box/self.nbins
                zposes.append(zpos)
            contig.append(zposes)
        return contig

    def find_nearest_halo(self):
        """Find the single most massive halos associated with absorption near a sightline, possibly via a subhalo."""
        (halos, subhalos) = self.find_nearby_halos()
        outhalos = np.zeros(self.NumLos,dtype=int)-1
        for ii in xrange(self.NumLos):
            subhalo_parent = list(self.sub_sub_index[subhalos[ii]])
            both = list(set(subhalo_parent+halos[ii]))
            if len(both) > 0:
                vir_vel = self.virial_vel(both)
                ind = np.where(vir_vel == np.max(vir_vel))
                outhalos[ii] = both[ind[0][0]]
            else:
                outhalos[ii] = -1
        return (outhalos, 0)

    def find_nearby_halos(self):
        """Find halos and subhalos associated with absorption near a sightline"""
        try:
            return (self.spectra_halos, self.spectra_subhalos)
        except AttributeError:
            pass
        zpos = self.get_contiguous_regions(thresh = 1e19, relthresh = 1e-2)
        (halos, _) = self.assign_to_halo(zpos, self.sub_radii, self.sub_cofm)
        (subhalos, _) = self.assign_to_halo(zpos, self.sub_sub_radii, self.sub_sub_cofm)
        #Merge absorption features inside the same halo
        for ii in xrange(self.NumLos):
            halos[ii] = list(set(halos[ii]))
            subhalos[ii] = list(set(subhalos[ii]))
        print("no. halos: ",sum([len(hh) for hh in halos])," mult halos: ",sum([len(hh) > 1 for hh in halos]))
        print("no. subhalos: ",sum([len(hh) for hh in subhalos])," mult subhalos: ",sum([len(hh) > 1 for hh in subhalos]))
        self.spectra_halos = halos
        self.spectra_subhalos = subhalos
        return (halos, subhalos)

    def get_stellar_mass_function(self):
        """Plot the galaxy stellar mass function for a snapshot."""
        subs=subfindhdf.SubFindHDF5(self.base, self.num)
        stellar_mass = subs.get_grp("GroupMassType")[:,4]*1e10/0.7
        #Could also use subhalo stellar mass: they are similar
        #stellar_mass = subs.get_sub("SubhaloMassType")[:,4]
        bins = np.logspace(6,12)
        dlogM = np.diff(np.log10(bins[:2]))
        volume = (25/0.7)**3
        (gsmf,sm) = np.histogram(stellar_mass, bins=bins)
        sm = (sm[1:]+sm[:-1])/2.
        gsmf = gsmf/volume/dlogM
        return (sm, gsmf)

    def _plot_breakdown(self, array, filt, low, high, labels, dv, log=True):
        """
        Helper function to plot something broken down by halo mass
        """
        #Find virial velocity
        (halo, _) = self.find_nearest_halo()
        ind = np.where(halo[filt] > 0)
        virial = self.virial_vel(halo[filt][ind])
        array = array[filt]
        #Make bins
        if log:
            func = plt.semilogx
            v_table = 10**np.arange(np.min(np.log10(array)),np.max(np.log10(array)) , dv)
        else:
            func = plt.plot
            v_table = np.arange(np.min(array),np.max(array) , dv)
        vbin = np.array([(v_table[i]+v_table[i+1])/2. for i in range(0,np.size(v_table)-1)])
        #Histogram of vel width
        vhist = np.histogram(array, v_table)[0]
        vhist[np.where(vhist == 0)] = 1
        colors = ("red", "purple", "cyan")
        lss = ("--", ":", "-")
        #Histogram of vel width for all halos in given virial velocity bin
        for ii in xrange(len(low)):
            vind = np.where((virial > low[ii])*(virial < high[ii]))
            vhist2 = np.histogram(array[ind][vind], v_table)[0]
            func(vbin, vhist2/(1.*vhist), color=colors[ii], ls=lss[ii], label=labels[ii])
#         vind = np.where(halo[filt] < 0)
#         vhist2 = np.histogram(array[vind], v_table)[0]
#         func(vbin, vhist2/(1.*vhist), color="grey", ls="-.", label="Field")

    def plot_Z_vs_mass(self,color="blue", color2="darkblue"):
        """Plot the correlation between mass and metallicity, with a fit"""
        (halo, _) = self.find_nearest_halo()
        ind = np.where(halo > 0)
        met = self.get_metallicity()[ind]
        mind = np.where(met > 1e-4)
        halo = halo[ind]
        mass = self.sub_mass[halo]
        mass = mass[mind]
        met = met[mind]
        self._plot_2d_contour(mass+0.1, met, 10, "Z mass", color, color2)
        plt.ylim(1e-4,1)

    def _plot_xx_vs_mass(self, xx, name = "xx", color="blue", color2="darkblue", log=True):
        """Helper function to plot something against virial velocity"""
        (halo, _) = self.find_nearest_halo()
        ind = np.where(halo > 0)
        halo = halo[ind]
        xx = xx[ind]
        virial = self.virial_vel(halo)+0.1
        self._plot_2d_contour(virial, xx, 10, name+" virial velocity", color, color2, ylog=log)


def combine_regions(condition, mindist=0):
    """Combine contiguous regions that are shorter than mindist"""
    reg = contiguous_regions(condition)
    #Find lengths to ignore
    if mindist > 0 and np.shape(reg)[0] > 1:
        newreg = np.array(reg[0,:])
        newreg.shape = (1,2)
        for ii in xrange(1,np.shape(reg)[0]):
            if reg[ii,0] - newreg[-1,1] < mindist:
                #Move the end point of the last segment to that of this one
                newreg[-1,1] = reg[ii,1]
            else:
                #This segment is far from the last one.
                #Add the new segment to the list
                newreg = np.vstack([newreg, reg[ii,:]])
        reg = newreg
    return reg

def contiguous_regions(condition):
    """Finds contiguous True regions of the boolean array "condition". Returns
    a 2D array where the first column is the start index of the region and the
    second column is the end index.
    If mindist != 0, ignores changes shorter than mindist
    """
    # Find the indicies of changes in "condition"
    d = np.diff(condition)
    idx, = d.nonzero()
    # We need to start things after the change in "condition". Therefore,
    # we'll shift the index by 1 to the right.
    idx += 1

    if condition[0]:
        # If the start of condition is True prepend a 0
        idx = np.r_[0, idx]

    if condition[-1]:
        # If the end of condition is True, append the length of the array
        idx = np.r_[idx, condition.size]

    # Reshape the result into two columns
    idx.shape = (-1,2)
    return idx
