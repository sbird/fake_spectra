"""Little module to find the path of a Cosmo box simulation"""

import os.path as path

base=path.expanduser("~/data/Cosmo/")

def get_name(sim, ff=True):
    """Get the directory for a simulation"""
    halo = "Cosmo"+str(sim)+"_V6"
    if ff:
        halo=path.join(halo,"L25n512")
    else:
        halo=path.join(halo,"L25n256")
    return path.join(base, halo)
