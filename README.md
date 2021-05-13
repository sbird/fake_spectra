= Flux extractor =

This is a small code for generating and analyzing simulated spectra from
Arepo/Gadget HDF5 simulation output. It is fast, parallel (Open-MP and MPI) and written in C++ and Python 3.
It really has two parts:
1) a C++/python 3 code which generates and analyses arbitrary spectra
2) A (slightly) maintained C++ command line program, extract, that generates
Lyman-alpha spectra and outputs them to a binary file.

If you are reading these instructions, I will assume you are using the first part:
the second program is for compatibility with older spectral extraction codes,
of which this is a rewrite, and so anyone who wants it should already know
how to use it.

== Installation ==

*This is a python 3 code*
*Make sure you install the python 3 libraries!*

Required Python libraries:
- numpy (core functionality)
- h5py (for saving)

Required C libraries:
- GSL

Optional libraries:
- matplotlib (if you want to plot)
- bigfile (to install, do 'pip install --user bigfile') for reading BigFile snapshot outputs from Yu Feng's MP-Gadget.
- mpi4py (Only if you want to use MPI feature in this package, it is very helpful for high resolution spectra in big simulations)

All these libraries can be installed with pip.

The easiest way to install the code is with pip:
```
pip3 install --user fake_spectra
```

If you need to use a pre-release version for any reason, you need this git repo.

First you need to check out the submodules:
```
git submodule update --init
```
Then compile it using:
```
python3 setup.py build
python3 setup.py install --user
```
On some systems you may have to add the directory it installs to
(usually $HOME/.local/lib) to your $PYTHONPATH

At time of writing, the code should compile with python2, once
python3-config in the Makefile is replaced with python2-config.
However, I do not guarantee bug-free operation, and strongly
recommend using python 3.

The test suite for the C++ module (only required during development)
requires Boost::Test and can be used with "make test"

== Usage ==
The main spectral generation routines are used can be called with:
```
import fake_spectra.spectra
spectra.Spectra(...)
```
Spectra takes two arguments: cofm and axis. If both are set to None,
the code will load them from a savefile. If there is no savefile you
must specify them. If N is the number of sightlines, axis has N entries.
axis specifies in which direction the sightline goes through the box.
Each entry of axis may be 1,2,3, with 1 being the x-axis, 2 being
the y-axis and 3 being the z-axis. In other words, axis is a
*1-indexed* array due to the original FORTRAN heritage of the code.

cofm is an Nx3 array. The column of cofm corresponding to axis is ignored:
it is left specified only so that you can easily generate sightlines
going through different axes in one call. Hence the following:

cofm = [[1,2,3], [4,5,6]], axis = 1

generates sightlines along the x-axis at y,z = 2,3 and y,z=5,6.
It produces identical output to

cofm = [[4,2,3], [1,5,6]], axis = 1

I have created a number of convenient wrappers for common configurations.
Two of these are:

randspectra.py - Generate spectra at random locations along the x axis.
Can optionally discard spectra that do not meet an HI column density threshold.
halospectra.py - Generate spectra through the center of halos.

Spectral generation routines take two arguments, base and num, which
specify where they should look for snapshot output. They will search:
`$(base)/snapdir_$(num)/snap_$(num).hdf5`
Note that num is padded with zeros to three characters, so passing '40' will result in '040'.

Column densities can be generated for arbitrary ions with the method get_col_density(elem, ion)
For neutral gas, pass ion 1. For the sum of all ionic species, pass ion -1.

Optical depths can be generated for arbitrary lines with the method get_tau(elem, ion, line)
Line data is loaded from a copy of atom.dat taken from VPFIT.

Thus, to generate randomly positioned Lyman-series spectra and associated HI column densities,
one would use this script:

```
from fake_spectra.randspectra import RandSpectra
# Only if you want the MPI feature; otherwise, MPI is None by default
from mpi4py import MPI

rr = RandSpectra(5, "MySim", MPI=MPI, thresh=0.)
rr.get_tau("H",1,1215)
#Lyman-beta
rr.get_tau("H",1,1025)
rr.get_col_density("H",1)
#Save spectra to file
rr.save_file()
```
If you want a set of spectra on a regular grid (`60*60` in this example), instead of `RandSpectra()` call :

```
from fake_spectra.griddedspectra import GriddedSpectra
# Only if you want the MPI feature; otherwise, MPI is None by default
from mpi4py import MPI

gs = GriddedSpectra(5, "MySim", nspec=60, MPI=MPI, thresh=0.0)
gs.get_tau("H",1,1025)
...

```

Note that the wavelength of the transition always rounds down,
so lyman alpha at 1215.67 A is 1215, not 1216!

In order to use the `MPI feature`, you need to run the code with `mpirun` or `mpiexec` or ... :

```
mpirun -n 10 python script_above.py

```

The MPI feature doese NOT work on Voronoi meshes due to dependency of particles and
the whole snapshot needs to be loaded on memory. So, for large simulations like Arepo 
to use the MPI feature you have to manually pass the `kernel="tophat"` to `RandSpectra()` 
or `GriddedSpectra()`. This choice of kernel is accurate. For more see [here](https://github.com/sbird/fake_spectra/blob/d0ec14e049fdced086fd9de9b9f16afddc1c71bc/fake_spectra/spectra.py#L76)

Generated spectra will be saved into HDF5 files, for ease of later analysis.
Each spectral generation routine saves spectra to a differently named file.

To load them again, use the PlottingSpectra routines:
```
from fake_spectra.plot_spectra import PlottingSpectra

ps = PlottingSpectra(num=5,base="MySim", label="My label",savefile="mysavefile.hdf5")
ps.plot_cddf("H",1)
```

You can also compute the temperature-density relation and mean IGM temperature with the tempdens module:

```
from fake_spectra.tempdens import fit_td_rel_plot

fit_td_rel_plot(5, "MySim", plot=True)
```
