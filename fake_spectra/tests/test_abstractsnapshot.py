# -*- coding: utf-8 -*-
"""Tests for the snapshot reading abstraction."""

import os
import h5py
import numpy as np
from fake_spectra import abstractsnapshot as absn

def test_hdf5_get_npart(tmp_path):
    """Check the particle count, including the high word, is read as an int64 array.
    Header attributes are uint32, as written by Gadget/Arepo."""
    snapdir = tmp_path / "snapdir_005"
    os.mkdir(snapdir)
    npart = np.array([10, 20, 0, 0, 5, 1], dtype=np.uint32)
    highword = np.array([1, 0, 0, 0, 2, 0], dtype=np.uint32)
    with h5py.File(snapdir / "snap_005.0.hdf5", "w") as ff:
        ff.create_group("Header")
        ff["Header"].attrs["NumPart_Total"] = npart
        ff["Header"].attrs["NumPart_Total_HighWord"] = highword
    snap = absn.HDF5Snapshot(5, str(tmp_path), None)
    got = snap.get_npart()
    expected = np.array([2**32 + 10, 20, 0, 0, 2 * 2**32 + 5, 1], dtype=np.int64)
    assert np.array_equal(got, expected)
