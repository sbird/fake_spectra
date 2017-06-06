"""python-config when used with homebrew on the mac
does not report the correct linker path.
Instead we have to open python and use distutils to get it."""

from __future__ import print_function
import distutils
import distutils.sysconfig

print(distutils.sysconfig.get_python_lib(standard_lib=True))
