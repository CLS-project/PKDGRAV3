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
import itertools


def align_structures(
    reference,
    generated,
    *,
    pos_field="rPot",
    dist2_thresh=1e-6,
    max_block=4,
    verbose=False,
):
    """
    Align `generated` to `reference` by:
      1) sorting both by descending number of DM, ascending rMax
      2) computing squared distances of aligned positions
      3) fixing local mis-orderings (small contiguous blocks)
      4) asserting remaining distances are small

    Returns:
        reference_sorted, generated_aligned
    """

    assert reference.dtype == generated.dtype
    assert reference.shape == generated.shape

    # ---- step 1: sort both by mass (descending)
    order_s = np.lexsort((
        reference["rMax"],
        -reference["nDM"],
    ))
    order_p = np.lexsort((
        generated["rMax"],
        -generated["nDM"],
    ))

    reference_s   = reference[order_s]
    generated_s = generated[order_p].copy()

    # ---- helpers
    def dist2(a, b):
        d = a[pos_field] - b[pos_field]
        return np.einsum("ij,ij->i", d, d)

    def contiguous_blocks(mask):
        idx = np.where(mask)[0]
        if len(idx) == 0:
            return []

        blocks = []
        start = prev = idx[0]
        for i in idx[1:]:
            if i == prev + 1:
                prev = i
            else:
                blocks.append((start, prev + 1))
                start = prev = i
        blocks.append((start, prev + 1))
        return blocks

    # ---- step 2: detect bad alignments
    d2 = dist2(reference_s, generated_s)
    bad = d2 > dist2_thresh

    if not np.any(bad):
        return reference_s, generated_s

    # ---- step 3: fix small local blocks
    blocks = contiguous_blocks(bad)

    for i, j in blocks:
        size = j - i
        if size > max_block:
            raise AssertionError(
                f"Alignment block too large ({size}) at [{i}:{j}); "
                "likely a real mismatch"
            )

        S = reference_s[i:j]
        P = generated_s[i:j]

        best_perm = None
        best_cost = np.inf

        for perm in itertools.permutations(range(size)):
            cost = 0.0
            for k, p in enumerate(perm):
                d = S[k][pos_field] - P[p][pos_field]
                cost += np.dot(d, d)
            if cost < best_cost:
                best_cost = cost
                best_perm = perm

        if verbose:
            print(f"Fixing block [{i}:{j}], best cost = {best_cost:.3e}")

        generated_s[i:j] = P[list(best_perm)]

    # ---- step 4: final validation
    d2_final = dist2(reference_s, generated_s)

    if np.any(d2_final > dist2_thresh):
        idx = np.argmax(d2_final)
        raise AssertionError(
            "Unresolved mismatch after alignment\n"
            f" index          : {idx}\n"
            f" reference rPot : {serial_s[idx][pos_field]}\n"
            f" generated rPot : {generated_s[idx][pos_field]}\n"
            f" dist2          : {d2_final[idx]:.3e}"
        )

    return reference_s, generated_s

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
@unittest.skipIf(not os.path.isfile('b0-final-ref.grp.npy'), "missing b0-final-ref.grp.npy")
@unittest.skipIf(not os.path.isfile('b0-final.fofstats.npy'), "missing b0-final.fofstats.npy")
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

        reference = np.fromfile("b0-final.fofstats.npy", dtype=TGT_dtype)
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = pathlib.Path(tmpdir)
            base = tmpdir / "fof.gs"
            msr.write_group_stats(f"{base}.serial", parallel=False)
            serial = np.fromfile(f"{base}.serial", dtype=TGT_dtype)
            # 
            # serial.tofile("b0-final.fofstats.npy")
            #
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

        # The centre of mass is more stable, but is relative to rPot
        reference['rcom'] += reference['rPot']
        serial['rcom'] += serial['rPot']

        reference,serial = align_structures(reference,serial)
        assert_structured_nearly_equal(reference,serial)

if __name__ == '__main__':
    print('Running test')
    unittest.main(verbosity=2,testRunner=xmlrunner.XMLTestRunner(output='test-reports'))
