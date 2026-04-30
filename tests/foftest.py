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

import numpy as np

def assert_structured_nearly_equal(a, b, *,
                                   rtol=1e-5,
                                   atol=1e-7,
                                   ignore_fields=("iGlobalGid", # assigned arbitrarily
                                                  "rPot",       # not stable
                                                  "rHalf",      # depends on center
                                                  "rMax",       # depends on center
                                                  "minPot"),    # not stable
                                   dtype_tols=None,):
    if dtype_tols is None:
        dtype_tols = {
            np.float32: dict(rtol=1e-5, atol=1e-7),
            np.float64: dict(rtol=1e-12, atol=1e-14),
        }

    assert a.dtype == b.dtype
    assert a.shape == b.shape

    for name in a.dtype.names:
        if name in ignore_fields:
            continue

        dt = a.dtype[name]

        # unwrap subarray dtype if present
        if dt.subdtype is not None:
            base_dtype, _ = dt.subdtype
        else:
            base_dtype = dt

        if np.issubdtype(base_dtype, np.floating):
            tol = dtype_tols.get(
                base_dtype.type,
                dict(rtol=1e-6, atol=1e-12)  # fallback
            )
            np.testing.assert_allclose(
                a[name], b[name],
                rtol=tol["rtol"], atol=tol["atol"],
                err_msg=f"Mismatch in field '{name}'"
            )
        else:
            np.testing.assert_array_equal(
                a[name], b[name],
                err_msg=f"Mismatch in field '{name}'"
            )

@ddt
@unittest.skipIf(not os.path.isfile('b0-final.std'), "missing b0-final.std")
@unittest.skipIf(not os.path.isfile('b0-final-ref.npz'), "missing b0-final-ref.npz")
class TestFof(unittest.TestCase):
    def testFof(self):
        # for the unit test these values are fixed!
        name = 'b0-final.std'
        dTau = 0.0004098
        nMinMembers = 10

        time = msr.load(name,bFindGroups = True,bMemGlobalGid = True,bPeriodic=True,bEwald=True,nReplicas=2,bEpsAccStep=True,bMemIntegerPosition=True,nMemEphemeral=8,nParaWrite=4)
        msr.domain_decompose()
        msr.build_tree(ewald=True)
        msr.fof(dTau,nMinMembers,statistics=False)
        msr.gravity(time=time,theta=0.5)
        msr.group_stats()

        b0ref = np.load("b0-final-ref.npz")
        reference = b0ref['fofstats']
        # reference = np.fromfile("b0-final-ref.fofstats.npz", dtype=TGT_dtype)
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
        a = b0ref['fof']
        self.assertTrue(np.all(np.where(a == 0)[0]==np.where(b == 0)[0]))
        (ua,ia,ra) = np.unique(a,return_index=True,return_inverse=True)
        m = b[ia]
        x = m[ra]

        np.testing.assert_array_equal(x, b)
        np.testing.assert_array_equal(serial['iGlobalGid'],unique_ids[1:])
        np.testing.assert_array_equal(serial['nDM'],counts[1:])

        # ua[i] -> m[i] = b[ia[i]]: old (reference grp) group ID to new (generated) group ID
        gid_map = dict(zip(ua[1:].tolist(), m[1:].tolist()))
        # new (generated) group ID -> index in serial/parallel stats
        serial_id_to_pos   = {int(gid): i for i, gid in enumerate(serial['iGlobalGid'])}
        parallel_id_to_pos = {int(gid): i for i, gid in enumerate(parallel['iGlobalGid'])}
        serial_perm   = np.array([serial_id_to_pos[gid_map[int(gid)]]   for gid in reference['iGlobalGid']])
        parallel_perm = np.array([parallel_id_to_pos[gid_map[int(gid)]] for gid in reference['iGlobalGid']])

        # Convert rcom from rPot-relative to absolute before comparing
        reference['rcom'] += reference['rPot']
        serial['rcom']    += serial['rPot']
        parallel['rcom']  += parallel['rPot']

        assert_structured_nearly_equal(reference, serial[serial_perm])
        assert_structured_nearly_equal(reference, parallel[parallel_perm])

if __name__ == '__main__':
    print('Running test')
    unittest.main(verbosity=2,testRunner=xmlrunner.XMLTestRunner(output='test-reports'))
