from cosmology import ClassCosmology, SimpleCosmology
import math
import unittest
from ddt import ddt, data, file_data, unpack
import xmlrunner
import os

@ddt
class TestCosmology(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.csm = SimpleCosmology(omega_matter=0.32)

    @data(1.0,0.1,0.01)
    def testTimeIdentity(self,a):
        self.assertAlmostEqual(self.csm.time_to_expansion(self.csm.expansion_to_time(a)),a,places=15)

    @data([1.0,0.32710340943409555],[0.1,0.012871250907720502],[0.01,0.00040716861571363553])
    @unpack
    def testExp2Time(self,a,t):
        self.assertAlmostEqual(self.csm.expansion_to_time(a),t,places=15)

    @data([1.0,0.32710340943409555],[0.1,0.012871250907720502],[0.01,0.00040716861571363553])
    @unpack
    def testTime2Exp(self,a,t):
        self.assertAlmostEqual(self.csm.time_to_expansion(t),a,places=15)

    def testExp2Hub(self):
        self.assertAlmostEqual(self.csm.expansion_to_hubble(1.0), math.sqrt(math.pi*8/3),places=15)
        self.assertAlmostEqual(self.csm.expansion_to_hubble(0.1),51.83167454117027,places=13)
        self.assertAlmostEqual(self.csm.expansion_to_hubble(0.01),1637.3244723688608,places=11)

    def testTime2Hub(self):
        self.assertAlmostEqual(self.csm.time_to_hubble(0.32710340943409555),math.sqrt(math.pi*8/3),places=15)
        self.assertAlmostEqual(self.csm.time_to_hubble(0.012871250907720502),51.83167454117027,places=13)
        self.assertAlmostEqual(self.csm.time_to_hubble(0.00040716861571363553),1637.3244723688608,places=11)

    def testComoveDriftFac(self):
        self.assertAlmostEqual(self.csm.comoving_drift_factor(0.1,0.01),0.05912002238902069,places=15)

    def testComoveKickFac(self):
        self.assertAlmostEqual(self.csm.comoving_kick_factor(0.1,0.01),0.024310129641862636,places=15)

    def testComoveGrowth(self):
        (D1,D2,D3a,D3b,D3c,f1,f2,f3a,f3b,f3c) = self.csm.comoving_growth(1.0)
        self.assertAlmostEqual(D1,  0.7906690062585626,places=10)
        self.assertAlmostEqual(D2,  0.2701133927580027,places=10)
        self.assertAlmostEqual(D3a, 0.16761178269550583,places=10)
        self.assertAlmostEqual(D3b, 0.23958097076482035,places=10)
        self.assertAlmostEqual(D3c, 0.07178807746250378,places=10)
        self.assertAlmostEqual(f1,  0.5318013708768904,places=10)
        self.assertAlmostEqual(f2,  1.0780894851672784,places=10)
        self.assertAlmostEqual(f3a, 1.6263116549213064,places=10)
        self.assertAlmostEqual(f3b, 1.6272836282587904,places=10)
        self.assertAlmostEqual(f3c, 1.625213026413398,places=10)
        (D1,D2,D3a,D3b,D3c,f1,f2,f3a,f3b,f3c) = self.csm.comoving_growth(0.1)
        self.assertAlmostEqual(D1,  0.09996140149711615,places=10)
        self.assertAlmostEqual(D2,  0.004282470274471333,places=10)
        self.assertAlmostEqual(D3a, 0.00033295780671554245,places=10)
        self.assertAlmostEqual(D3b, 0.00047565452931935787,places=10)
        self.assertAlmostEqual(D3c, 0.00014269602995431666,places=10)
        self.assertAlmostEqual(f1,  0.9988427590347082,places=10)
        self.assertAlmostEqual(f2,  1.9977298409492161,places=10)
        self.assertAlmostEqual(f3a, 2.9966204670359717,places=10)
        self.assertAlmostEqual(f3b, 2.99662356330118,places=10)
        self.assertAlmostEqual(f3c, 2.9966170570829243,places=10)

# This needs an HDF5 file
@unittest.skipIf(not os.path.isfile('euclid_flagship_500.hdf5'), "missing euclid_flagship_500.hdf5")
class TestCosmologyClass(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.csm1 = ClassCosmology('euclid_flagship_500.hdf5',L=1,As=2.1e-09,ns=0.96)
        self.csm2 = SimpleCosmology(omega_matter=self.csm1.omega_matter + (0.0011978 + 0.000208644 ),
                                     omega_lambda=self.csm1.omega_lambda,
                                     omega_radiation=self.csm1.omega_radiation + 1.27031e-05,
                                     omega_dark_energy=self.csm1.omega_dark_energy,
                                     w0=self.csm1.w0,wa=self.csm1.wa)

    def testHubble(self):
        self.assertAlmostEqual(self.csm1.expansion_to_hubble(a=1.0), self.csm2.expansion_to_hubble(a=1.0),places=15)
        self.assertAlmostEqual(self.csm1.expansion_to_hubble(a=0.1), self.csm2.expansion_to_hubble(a=0.1),delta=0.005)

if __name__ == '__main__':
    print('Running test')
    unittest.main(verbosity=2,testRunner=xmlrunner.XMLTestRunner(output='test-reports'))
