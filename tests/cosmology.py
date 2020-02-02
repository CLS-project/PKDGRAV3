from __future__ import division
import MASTER as MSR
from CSM import CSM
import math
import unittest
from ddt import ddt, data, file_data, unpack
import xmlrunner

@ddt
class TestCosmology(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.csm = CSM()
        cls.csm.SetCosmology(dHubble0=math.sqrt(math.pi*8/3),dOmega0=0.32,dSigma8=0.83,ns=0.96)

    @data(1.0,0.1,0.01)
    def testTimeIdentity(self,a):
        self.assertAlmostEqual(self.csm.Time2Exp(self.csm.Exp2Time(a)),a,places=15)

    @data([1.0,0.27718310726907297],[0.1,0.012137233682214097],[0.01,0.0004046022074540854])
    @unpack
    def testExp2Time(self,a,t):
        self.assertAlmostEqual(self.csm.Exp2Time(a),t,places=15)

    @data([1.0,0.27718310726907297],[0.1,0.012137233682214097],[0.01,0.0004046022074540854])
    @unpack
    def testTime2Exp(self,a,t):
        self.assertAlmostEqual(self.csm.Time2Exp(t),a,places=15)

    def testExp2Hub(self):
        self.assertAlmostEqual(self.csm.Exp2Hub(1.0), math.sqrt(math.pi*8/3),places=15)
        self.assertAlmostEqual(self.csm.Exp2Hub(0.1),57.0131668907651,places=13)
        self.assertAlmostEqual(self.csm.Exp2Hub(0.01),1654.62783665947,places=11)

    def testTime2Hub(self):
        self.assertAlmostEqual(self.csm.Time2Hub(0.27718310726907297),math.sqrt(math.pi*8/3),places=15)
        self.assertAlmostEqual(self.csm.Time2Hub(0.012137233682214097),57.0131668907651,places=13)
        self.assertAlmostEqual(self.csm.Time2Hub(0.0004046022074540854),1654.62783665948,places=11)

    def testComoveDriftFac(self):
        self.assertAlmostEqual(self.csm.ComoveDriftFac(0.1,0.01),0.04525620448976383,places=15)

    def testComoveKickFac(self):
        self.assertAlmostEqual(self.csm.ComoveKickFac(0.1,0.01),0.021268935133646454,places=15)

    def testComoveGrowth(self):
        (D1,D2,f1,f2) = self.csm.ComoveGrowth(1.0)
        self.assertAlmostEqual(D1, 0.47863590415639595,places=15)
        self.assertAlmostEqual(D2,-0.10199487849270733,places=15)
        self.assertAlmostEqual(f1, 0.5114165755213692,places=15)
        self.assertAlmostEqual(f2, 1.0461127333263733,places=15)
        (D1,D2,f1,f2) = self.csm.ComoveGrowth(0.1)
        self.assertAlmostEqual(D1, 0.08934869147265322,places=15)
        self.assertAlmostEqual(D2,-0.0034426240751702234,places=15)
        self.assertAlmostEqual(f1, 0.8952794949702897,places=15)
        self.assertAlmostEqual(f2, 1.7962589373634128,places=15)

# # This needs an HDF5 file
# class TestCosmologyClass(unittest.TestCase):
#     @classmethod
#     def setUpClass(self):
#         self.csm1 = CSM()
#         self.csm1.ClassRead('euclid_flagship_500.hdf5',Lbox=1,As=2.1e-09,ns=0.96)
#         self.csm2 = CSM()
#         self.csm2.SetCosmology(dHubble0=math.sqrt(math.pi*8/3),As=2.1e-09,ns=0.96,
#         		dOmega0=self.csm1.dOmega0 + (0.0011978 + 0.000208644 ),
#         		dLambda=self.csm1.dLambda,
#         		dOmegaRad=self.csm1.dOmegaRad + 1.27031e-05,
#         		dOmegaDE=self.csm1.dOmegaDE,w0=self.csm1.w0,wa=self.csm1.wa)

#     def testHubble(self):
#         self.assertAlmostEqual(self.csm1.Exp2Hub(a=1.0), self.csm2.Exp2Hub(a=1.0),places=15)
#         self.assertAlmostEqual(self.csm1.Exp2Hub(a=0.1), self.csm2.Exp2Hub(a=0.1),delta=0.005)

if __name__ == '__main__':
    print('Running test')
    unittest.main(verbosity=2,testRunner=xmlrunner.XMLTestRunner(output='test-reports'))
