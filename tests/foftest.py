from __future__ import division
import os
from math import sqrt,pi
import numpy as np
import unittest
from ddt import ddt, data, file_data, unpack
import xmlrunner
from MASTER import MSR # Once MSR is imported, simulation mode is no longer automatically entered

@ddt
@unittest.skipIf(not os.path.isfile('b0-final.std'), "missing b0-final.std")
@unittest.skipIf(not os.path.isfile('b0-final-ref.grp.npy'), "missing b0-final-ref.grp.npy")
class TestFof(unittest.TestCase):
    def testFof(self):
        # for the unit test these values are fixed!
        name = 'b0-final.std'
        dTau = 0.0004098
        nMinMembers = 10
        msr = MSR()
        msr.setParameters(bFindGroups = True,bMemGlobalGid = True)
        time = msr.Load(name)
        msr.DomainDecomp()
        msr.BuildTree()
        msr.Fof(dTau,nMinMembers)
        msr.Reorder()
        b = msr.GetArray(field=16,time=time) # 16 = GLOBAL Group id (oGlobalGid)
        a = np.load('b0-final-ref.grp.npy')
        (ua,ia,ra) = np.unique(a,return_index=True,return_inverse=True)
        m = b[ia]
        x = m[ra]
        self.assertTrue(np.array_equal(x,b))

if __name__ == '__main__':
    print('Running test')
    unittest.main(verbosity=2,testRunner=xmlrunner.XMLTestRunner(output='test-reports'))
