from __future__ import division
import MASTER as MSR
from CSM import CSM
import math
import unittest
import xmlrunner

class TestCosmology(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.csm = CSM()
        self.csm.SetCosmology(dHubble0=math.sqrt(math.pi*8/3),dOmega0=0.32,dSigma8=0.83,ns=0.96)

    def testHubble(self):
        self.assertEqual(self.csm.Exp2Hub(1.0), math.sqrt(math.pi*8/3))
        self.assertEqual(self.csm.Exp2Hub(0.1),57.013166890765135)
        self.assertEqual(self.csm.Exp2Hub(0.01),1654.627836659466)

    def testExp2Time(self):
        self.assertEqual(self.csm.Exp2Time(1.0),0.27718310726907297)
        self.assertEqual(self.csm.Exp2Time(0.1),0.012137233682214097)
        self.assertEqual(self.csm.Exp2Time(0.01),0.0004046022074540854)

    def testTime2Exp(self):
        self.assertAlmostEqual(self.csm.Time2Exp(0.27718310726907297),1.0,places=15)
        self.assertAlmostEqual(self.csm.Time2Exp(0.012137233682214097),0.1,places=15)
        self.assertAlmostEqual(self.csm.Time2Exp(0.0004046022074540854),0.01,places=15)

if __name__ == '__main__':
    print('Running test')
    unittest.main(verbosity=2,testRunner=xmlrunner.XMLTestRunner(output='test-reports'))
