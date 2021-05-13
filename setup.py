from distutils.core import setup
from distutils.extension import Extension
import distutils.sysconfig
import numpy
import tempfile
import os
import subprocess
import shutil

def check_for_openmp():
    """Check  whether the default compiler supports OpenMP.

    This routine is adapted from yt, thanks to Nathan
    Goldbaum. See https://github.com/pynbody/pynbody/issues/124"""

    if os.getenv('DISABLE_OPENMP') is not None:
        return False

    # Create a temporary directory
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    # Get compiler invocation
    compiler = os.environ.get('CC',
                              distutils.sysconfig.get_config_var('CC'))

    # make sure to use just the compiler name without flags
    compiler = compiler.split()[0]

    # Attempt to compile a test script.
    # See http://openmp.org/wp/openmp-compilers/
    filename = r'test.c'
    with open(filename,'w') as f :
        f.write(
        "#include <omp.h>\n"
        "#include <stdio.h>\n"
        "int main() {\n"
        "#pragma omp parallel\n"
        "printf(\"Hello from thread %d, nthreads %d\\n\", omp_get_thread_num(), omp_get_num_threads());\n"
        "}"
        )

    try:
        with open(os.devnull, 'w') as fnull:
            exit_code = subprocess.call([compiler, '-fopenmp', filename],
                                        stdout=fnull, stderr=fnull)
    except OSError :
        exit_code = 1

    # Clean up
    os.chdir(curdir)
    shutil.rmtree(tmpdir)

    if exit_code == 0:
        return True
    else:
        import multiprocessing, platform
        cpus = multiprocessing.cpu_count()
        if cpus>1:
            print ("""WARNING

OpenMP support is not available in your default C compiler, even though
your machine has more than one core available.

Some routines in pynbody are parallelized using OpenMP and these will
only run on one core with your current configuration.
""")
            if platform.uname()[0]=='Darwin':
                print ("""Since you are running on Mac OS, it's likely that the problem here
is Apple's Clang, which does not support OpenMP at all. The easiest
way to get around this is to download the latest version of gcc from
here: http://hpc.sourceforge.net. After downloading, just point the
CC environment variable to the real gcc and OpenMP support should
get enabled automatically. Something like this -

sudo tar -xzf /path/to/download.tar.gz /
export CC='/usr/local/bin/gcc'
python setup.py clean
python setup.py build

""")
            print ("""Continuing your build without OpenMP...\n""")

        return False

extra_compile_args=['-ffast-math',]
extra_link_args = ['-ffast-math',]

try:
    gsl_libs = subprocess.check_output(["gsl-config", "--libs"], stderr=subprocess.STDOUT, universal_newlines=True)
    extra_link_args += gsl_libs.split()
    gsl_incl = subprocess.check_output(["gsl-config", "--cflags"], stderr=subprocess.STDOUT, universal_newlines=True)
    extra_compile_args += gsl_incl.split()
except subprocess.CalledProcessError as e:
    print(e.output)
    raise

if check_for_openmp():
    extra_compile_args += ['-fopenmp',]
    #gcc specific
    extra_link_args += ['-lgomp',]

cmodule = [
        Extension("fake_spectra._spectra_priv",
            ["fake_spectra/py_module.cpp",
             "fake_spectra/absorption.cpp",
             "fake_spectra/index_table.cpp",
             "fake_spectra/Faddeeva.cpp",
             "fake_spectra/part_int.cpp",
            ],
            depends = [
                "fake_spectra/Faddeeva.h",
                "fake_spectra/absorption.h",
                "fake_spectra/index_table.h",
                "fake_spectra/part_int.h",
                "fake_spectra/singleabs.h",]
            ,
            extra_compile_args=extra_compile_args,
            extra_link_args=extra_link_args,
            include_dirs = ["fake_spectra/", numpy.get_include()])]

setup(
    name="fake_spectra",
    version='2.2.0',
    author="Simeon Bird",
    author_email="spb@ucr.edu",
    #Use the subclass which adds openmp flags as appropriate
#     cmdclass = {'build_ext': build_ext_subclass },
    url="http://github.com/sbird/fake_spectra",
    description="Analysis tools for generating artificial spectra from simulations.",
    packages = ['fake_spectra', 'fake_spectra.tests', 'fake_spectra.cloudy_tables'],
    requires=['numpy', 'h5py','scipy'],
    package_data = {
            'fake_spectra.tests': ['*.npz'],
            'fake_spectra': ['data/TREECOOL*', '*.dat'],
            'fake_spectra.cloudy_tables': ['ion_out_*/cloudy_table.npz']
           },
    ext_modules = cmodule,
    classifiers = ["Development Status :: 4 - Beta",
                   "Intended Audience :: Developers",
                   "Intended Audience :: Science/Research",
                   "License :: OSI Approved :: MIT License",
                   "Programming Language :: Python :: 3",
                   "Topic :: Scientific/Engineering :: Astronomy",
                   "Topic :: Scientific/Engineering :: Visualization"]
)
