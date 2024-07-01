from __future__ import division
import os,sys
from math import sqrt,pi
import numpy as np
import unittest
from ddt import ddt, data, file_data, unpack
import xmlrunner
import PKDGRAV as msr # Once PKDGRAV is imported, simulation mode is no longer automatically entered

@ddt
@unittest.skipIf(not os.path.isfile('b0-final.std'), "missing b0-final.std")
@unittest.skipIf(not os.path.isfile('b0-final-ref.grp.npy'), "missing b0-final-ref.grp.npy")
class TestFof(unittest.TestCase):
    def testFof(self):
        # for the unit test these values are fixed!
        name = 'b0-final.std'
        dTau = 0.0004098
        nMinMembers = 10
        time = msr.load(name,bFindGroups = True,bMemGlobalGid = True,nMemEphemeral=8)
        msr.domain_decompose()
        msr.build_tree()
        msr.fof(dTau,nMinMembers)
        msr.reorder()
        b = msr.get_array(field=msr.FIELD_GLOBAL_GID,time=time) # 16 = GLOBAL Group id (oGlobalGid)
        a = np.load('b0-final-ref.grp.npy')
        self.assertTrue(np.all(np.where(a == 0)[0]==np.where(b == 0)[0]))
        (ua,ia,ra) = np.unique(a,return_index=True,return_inverse=True)
        m = b[ia]
        x = m[ra]
        self.assertTrue(np.array_equal(x,b))

if __name__ == '__main__':
    print('Running test')
    unittest.main(verbosity=2,testRunner=xmlrunner.XMLTestRunner(output='test-reports'))
