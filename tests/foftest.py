from __future__ import division
import os,sys
from math import sqrt,pi
import tempfile
import pathlib
import numpy as np
import unittest
from ddt import ddt, data, file_data, unpack
import xmlrunner
import PKDGRAV as msr # Once PKDGRAV is imported, simulation mode is no longer automatically entered
from PKDGRAV import PKD_FIELD

TGT_dtype = np.dtype([
    ("rPot",               np.float32, (3,)),
    ("minPot",             np.float32),
    ("rcen",               np.float32, (3,)),
    ("rcom",               np.float32, (3,)),
    ("vcom",               np.float32, (3,)),
    ("angular",            np.float32, (3,)),
    ("inertia",            np.float32, (6,)),
    ("sigma",              np.float32),
    ("rMax",               np.float32),
    ("fMass",              np.float32),
    ("fEnvironDensity0",   np.float32),
    ("fEnvironDensity1",   np.float32),
    ("rHalf",              np.float32),
    ("nBH",                np.int32),
    ("nStar",              np.int32),
    ("nGas",               np.int32),
    ("nDM",                np.int32),
    ("iGlobalGid",         np.uint64),
], align=False)

@ddt
@unittest.skipIf(not os.path.isfile('b0-final.std'), "missing b0-final.std")
@unittest.skipIf(not os.path.isfile('b0-final-ref.grp.npy'), "missing b0-final-ref.grp.npy")
class TestFof(unittest.TestCase):
    def testFof(self):
        # for the unit test these values are fixed!
        name = 'b0-final.std'
        dTau = 0.0004098
        nMinMembers = 10
        time = msr.load(name,bFindGroups = True,bMemGlobalGid = True,nMemEphemeral=8,nParaWrite=4)
        msr.domain_decompose()
        msr.build_tree()
        msr.fof(dTau,nMinMembers)

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = pathlib.Path(tmpdir)
            base = tmpdir / "fof.gs"
            msr.write_group_stats(f"{base}.serial", parallel=False)
            serial = np.fromfile(f"{base}.serial", dtype=TGT_dtype)
            msr.write_group_stats(f"{base}.parallel", parallel=4)
            parallel = np.concatenate([
                np.fromfile(f"{base}.parallel.{i}", dtype=TGT_dtype)
                for i in range(4)
            ])

        msr.reorder()
        b = msr.get_array(field=PKD_FIELD.FIELD_GLOBAL_GID,time=time) # 16 = GLOBAL Group id (oGlobalGid)
        unique_ids, counts = np.unique(b, return_counts=True)
        a = np.load('b0-final-ref.grp.npy')
        self.assertTrue(np.all(np.where(a == 0)[0]==np.where(b == 0)[0]))
        (ua,ia,ra) = np.unique(a,return_index=True,return_inverse=True)
        m = b[ia]
        x = m[ra]
        np.testing.assert_array_equal(x, b)
        np.testing.assert_array_equal(serial, parallel)
        np.testing.assert_array_equal(serial['iGlobalGid'],unique_ids[1:])
        np.testing.assert_array_equal(serial['nDM'],counts[1:])

if __name__ == '__main__':
    print('Running test')
    unittest.main(verbosity=2,testRunner=xmlrunner.XMLTestRunner(output='test-reports'))
