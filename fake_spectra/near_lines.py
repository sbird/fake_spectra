# -*- coding: utf-8 -*-
"""Find the particles which are near a sightline.

This replaces the C++ IndexTable used by the old _near_lines module.
The trick is that the caller only needs to know *whether* a particle is
within its smoothing length of some sightline, not which sightlines those
are. That makes it a nearest-neighbour query rather than a range query,
and scipy's cKDTree does nearest-neighbour queries in threaded C.

A KD-tree query for every particle is wasteful when (as is usual) almost
no particles are near a sightline, so particles are first passed through a
chaining mesh: a coarse boolean grid of cells which are within reach of a
sightline. Particles are binned by smoothing length in powers of two so
that a single mesh dilation radius covers each bin.
"""

import numpy as np
from scipy.spatial import cKDTree
from scipy.ndimage import maximum_filter

#The two coordinates perpendicular to each (1-indexed) sightline axis.
_PERP = {1: (1, 2), 2: (0, 2), 3: (0, 1)}

def near_lines(box, pos, hh, axis, cofm, ncell=1024, maxbin=5, workers=-1):
    """Find particles within a smoothing length of a sightline.

    The distance is the periodic distance in the two coordinates
    perpendicular to the sightline axis. Positions are assumed to be
    within [0, box), as they are in a snapshot.

    Arguments:
        box - periodic box size, in the same units as pos and cofm.
        pos - (npart, 3) particle positions.
        hh - (npart) particle smoothing lengths.
        axis - (nlos) sightline axes, 1-indexed, so between 1 and 3.
        cofm - (nlos, 3) sightline centres.
        ncell - side length of the chaining mesh. The default resolves
                a 25 Mpc box into 24 kpc cells. Larger is not better:
                the mesh lookup becomes memory bound well before it
                becomes more selective.
        maxbin - smoothing lengths larger than 2**maxbin cells skip the
                mesh and go straight to the tree. Without this a few
                large low-density particles force a very wide dilation.
        workers - threads for the tree query; -1 uses all cores.
    Returns:
        Sorted indices of the particles near a sightline.
    """
    npart = np.shape(pos)[0]
    near = np.zeros(npart, dtype=bool)
    if npart == 0 or np.size(axis) == 0:
        return np.nonzero(near)[0].astype(np.int32)
    cell = np.float32(box/ncell)
    #Bin the particles by smoothing length in powers of two.
    hbin = np.log2(np.maximum(hh/cell, np.float32(1.)))
    np.ceil(hbin, out=hbin)
    hbin = hbin.astype(np.int32)
    np.clip(hbin, 0, maxbin, out=hbin)
    nbin = int(hbin.max()) + 1
    #Bins with no particles in them need no mesh.
    used = np.zeros(nbin, dtype=bool)
    used[np.unique(hbin)] = True
    for ax in np.unique(axis):
        (aa, bb) = _PERP[int(ax)]
        lines = np.mod(cofm[axis == ax][:, [aa, bb]], box)
        #cKDTree wants its periodic points within [0, box).
        tree = cKDTree(lines, boxsize=box)
        mesh = _line_mesh(lines, cell, ncell, nbin, used, maxbin)
        ip = (pos[:, aa]*np.float32(1./cell)).astype(np.int32)
        np.clip(ip, 0, ncell-1, out=ip)
        iq = (pos[:, bb]*np.float32(1./cell)).astype(np.int32)
        np.clip(iq, 0, ncell-1, out=iq)
        maybe = mesh[hbin, ip, iq]
        del ip, iq, mesh
        #No need to look again at particles already known to be near a line.
        maybe &= np.logical_not(near)
        (cand,) = np.nonzero(maybe)
        if np.size(cand) == 0:
            continue
        perp = np.mod(pos[cand][:, [aa, bb]].astype(np.float64), box)
        #The nearest line is within h if and only if any line is.
        (dist, _) = tree.query(perp, k=1, workers=workers)
        near[cand[dist <= hh[cand]]] = True
    return np.nonzero(near)[0].astype(np.int32)

def _line_mesh(lines, cell, ncell, nbin, used, maxbin):
    """Build, for each smoothing length bin, a boolean mesh of the cells
    which are within reach of a sightline."""
    icell = np.clip((lines/cell).astype(np.intp), 0, ncell-1)
    occupied = np.zeros((ncell, ncell), dtype=bool)
    occupied[icell[:, 0], icell[:, 1]] = True
    mesh = np.empty((nbin, ncell, ncell), dtype=bool)
    for hb in range(nbin):
        if not used[hb]:
            mesh[hb] = False
            continue
        #A particle in this bin has h <= 2**hb cells, so a line more than
        #2**hb + 1 cells away in either coordinate cannot be within h.
        rr = 2**hb + 1
        if hb == maxbin or 2*rr + 1 >= ncell:
            mesh[hb] = True
        else:
            mesh[hb] = _dilate(icell, occupied, rr, ncell)
    return mesh

def _dilate(icell, occupied, rr, ncell):
    """Mark all mesh cells within rr cells of an occupied cell, periodically."""
    nlines = np.shape(icell)[0]
    #Scattering each line into its neighbourhood costs O(nlines rr^2), which
    #beats filtering the whole mesh only while the lines are sparse.
    if nlines*(2*rr + 1)**2 > 4*ncell*ncell:
        return maximum_filter(occupied, size=2*rr + 1, mode='wrap')
    grid = np.zeros((ncell, ncell), dtype=bool)
    off = np.arange(-rr, rr+1)
    #Chunked so the broadcast index arrays stay a sensible size.
    step = max(1, int(4e6 // (2*rr + 1)**2))
    for ss in range(0, nlines, step):
        chunk = icell[ss:ss+step]
        dx = np.mod(chunk[:, 0, None] + off, ncell)
        dy = np.mod(chunk[:, 1, None] + off, ncell)
        grid[dx[:, :, None], dy[:, None, :]] = True
    return grid
