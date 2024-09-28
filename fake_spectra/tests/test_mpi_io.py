"""
Unit tests for MPI I/O of the bigfile snapshots

These tests should pass on any bigfile snapshot (the pathis 
specifed with `base`) with any number of particles and any 
number of MPI ranks.
"""
import os
import unittest
import numpy as np
from mpi4py import MPI
from fake_spectra.abstractsnapshot import AbstractSnapshotFactory

class TestAbstractSnapshot(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.comm = MPI.COMM_WORLD
        cls.size = cls.comm.Get_size()
        cls.rank = cls.comm.Get_rank()

        cls.base = '../examle_bigfile/'
        cls.snap = 272
        cls.part_type = 0

        cls.abs_snap = AbstractSnapshotFactory(cls.snap, cls.base, MPI=MPI, log_level='debug')

    def test_bigfile_header(self):
        """Test the bigfile Header
        It is opened on rank 0 and broadcasted to all ranks.
        """
        header = self.abs_snap.header
        self.assertIsNotNone(header)
        self.assertIsInstance(header, dict)

    def test_npart(self):
        """Test the number of particles.
        It is retrieved from the bigfile Header."""
        npart = self.abs_snap.get_npart()
        self.assertIsNotNone(npart)
        print(f'rank = {self.rank} | TotNumPart = {npart}', flush=True)

    def test_block_path(self):
        """Test the block path retrieval."""
        block_path = self.abs_snap._get_block_path(self.part_type, 'Position')
        self.assertIsNotNone(block_path)
        print(f"rank = {self.rank} | block path = {block_path}", flush=True)

    def test_block_header(self):
        """Test the block header retrieval."""
        header = self.abs_snap._get_block_header_attr(self.part_type, 'Position')
        self.assertTrue( 'dtype' in header.keys())
        self.assertTrue( 'nmemb' in header.keys())
        self.assertEqual( len(header['blob_names']), 
                         len(header['blob_part_start']), 
                         len(header['blob_part_end']))

    def test_segmenting(self):
        """Test segmenting of the data."""
        nseg, chunk_size = self.abs_snap.get_n_segments(self.part_type)
        self.assertIsNotNone(nseg)
        self.assertIsNotNone(chunk_size)
        npart = self.abs_snap.get_npart()[self.part_type]
        # 1. Make sure the total number of particles distributed among all ranks
        # is equal to the total number of particles in the snapshot.
        self.assertEqual(
        np.sum(self.abs_snap.parts_rank), npart,
        f"sum(parts_rank) != npart[{self.part_type}] ({npart})"
        )
        # 2. make sure particle counts on all segments add up to the 
        # total number of particles on each rank
        self.assertGreaterEqual(
        nseg * chunk_size, self.abs_snap.parts_rank[self.rank],
        f"The last chunk should be smaller than or equal to the rest. Segmentaion on a single rank is wrong."
        )

        print(f'rank = {self.rank} | nseg: {nseg}, chunk_size: {chunk_size}', flush=True)
        print(f'rank = {self.rank} | particle laod : {self.abs_snap.parts_rank}', flush=True)

    def test_part_list(self):
        """Test particle list for each segment."""
        nseg, chunk_size = self.abs_snap.get_n_segments(self.part_type)
        segments = [0, 1, nseg - 1] # Don't change this
        start_seg, end_seg = [], []
        for i in segments:
            seg_start_end = self.abs_snap._segment_to_partlist(self.part_type, i)
            start_seg.append(seg_start_end[0])
            end_seg.append(seg_start_end[1])
        ## 3 tests:
        #1. check the continuity of the segments, the end is exclusive
        self.assertEqual(start_seg[1], end_seg[0])
        #2. check the size of the segments, i.e. <= chunk_size
        for i in range(len(segments)):
            self.assertLessEqual(end_seg[i] - start_seg[i], chunk_size)
        #3. check of the sum of segments is equal to the total particles on this rank
        self.assertEqual(
            end_seg[-1] - start_seg[0], self.abs_snap.parts_rank[self.rank],
            f"sum of segments != parts_rank[{self.rank}] ({self.abs_snap.parts_rank[self.rank]})"
        )

    def test_blobs_for_each_segment(self):
        """Test blob file names for a segment."""
        nseg, chunk_size = self.abs_snap.get_n_segments(self.part_type)
        segments = [0, 1, nseg-10, nseg-2, nseg - 1] # Don't change this
        for i in segments:
            # Type check
            self.assertTrue((start_blob.dtype == int) & (end_blob.dtype == int))

            seg_start, seg_end = self.abs_snap._segment_to_partlist(self.part_type, i)
            blobs, blob_paths, start_blob, end_blob = self.abs_snap._get_blobs_for_rank(
                segment=i, part_type=self.part_type, blockname='Position'
            )
            self.assertGreater(len(blobs), 0, "No blob files found.")
            for b in range(len(blobs)):
                # check the existence of the blob file
                self.assertTrue(os.path.exists(blob_paths[b]), f"blob file {blob_paths[b]} does not exist.")
                # The particle count in the blob file should be less than or equal to the segment size.
                self.assertLessEqual(
                    end_blob[b] - start_blob[b], seg_end - seg_start,
                    f"blob file {blob_paths[b]} has more particles than the segment."
                )

    def test_data(self):
        """Test data retrieval."""
        header = self.abs_snap._get_block_header_attr(self.part_type, 'Position')
        nseg, _ = self.abs_snap.get_n_segments(self.part_type)
        indices = [0, 8, nseg - 2, nseg-1]
        for i in indices:
            seg_start, seg_end = self.abs_snap._segment_to_partlist(self.part_type, i)
            data = self.abs_snap.get_data(self.part_type, 'Position', segment=i)
            self.assertIsNotNone(data)
            self.assertIsInstance(data, np.ndarray)

            #print(f' rank = {self.rank} | segment = {i} | data shape = {data.shape}  | {(seg_end - seg_start)}', flush=True)

            self.assertEqual(data.shape, (seg_end - seg_start, header['nmemb']))

if __name__ == '__main__':
    unittest.main()
