from __future__ import division
import os
from math import sqrt,pi
import numpy as np
import unittest
from ddt import ddt, data, file_data, unpack
import xmlrunner
from MASTER import MSR
from CSM import CSM

# @ddt
# class TestGravityBasic(unittest.TestCase):
#     @classmethod
#     def setUpClass(cls):
#         cls.csm = CSM()
#         cls.csm.ClassRead('euclid_flagship_500.hdf5',LinSpecies=['ncdm[0]','g','metric'],Lbox=1600,As=2.1e-09,ns=0.96)
#         cls.msr = MSR()
#         cls.time = cls.msr.GenerateIC(cls.csm,z=49,L=1600,grid=32,seed=43857)

#     def XtestGravityBasic(self):
#         self.msr.DomainDecomp()
#         self.msr.BuildTree()
#         self.msr.Gravity(time=self.time)

@ddt
@unittest.skipIf(not os.path.isfile('b0-final.std'), "missing b0-final.std")
@unittest.skipIf(not os.path.isfile('b0-final-np-asym-k1.acc.npy'), "missing b0-final-np-asym-k1.acc.npy")
class TestGravityB0Final(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.csm = CSM(dOmega0=0.32,dLambda=0.68,dSigma8=0.83,ns=0.96)
        cls.msr = MSR()
        cls.time = cls.msr.Load('b0-final.std')
        cls.a = np.load('b0-final-np-asym-k1.acc.npy')
        cls.maga = np.linalg.norm(cls.a,axis=1)

    @data([0.70,0.001],[0.60,0.0005],[0.55,0.0004],[0.40,8e-5],)
    @unpack
    def testGravityNonPeriodic(self,theta,target_rms):
        self.msr.setParameters(bPeriodic=False,bEwald=False,nReplicas=None)
        self.msr.DomainDecomp()
        self.msr.BuildTree()
        self.msr.Gravity(time=self.time,theta=theta) 
        self.msr.Reorder()
        a=self.msr.GetArray(field=2,time=self.time)
        relerr = np.linalg.norm(a-self.a,axis=1) / self.maga
        rms = np.std(relerr)
        print('rms',rms)
        self.assertLess(rms,target_rms)

@ddt
@unittest.skipIf(not os.path.isfile('b0-final.std'), "missing b0-final.std")
@unittest.skipIf(not os.path.isfile('b0-final-p0.10-asym-k1.acc.npy'), "missing b0-final-p0.10-asym-k1.acc.npy")
class TestGravityB0Periodic(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.csm = CSM(dOmega0=0.32,dLambda=0.68,dSigma8=0.83,ns=0.96)
        cls.msr = MSR()
        cls.time = cls.msr.Load('b0-final.std')
        cls.a = np.load('b0-final-p0.10-asym-k1.acc.npy')
        cls.maga = np.linalg.norm(cls.a,axis=1)

    @data([0.70,0.0011],[0.60,0.0005],[0.55,0.0004],[0.40,8.1e-5],)
    @unpack
    def testGravityPeriodic(self,theta,target_rms):
        self.msr.setParameters(bPeriodic=True,bEwald=True,nReplicas=2)
        self.msr.DomainDecomp()
        self.msr.BuildTree()
        self.msr.Gravity(time=self.time,theta=theta) 
        self.msr.Reorder()
        a=self.msr.GetArray(field=2,time=self.time)
        relerr = np.linalg.norm(a-self.a,axis=1) / self.maga
        rms = np.std(relerr)
        print('rms',rms)
        self.assertLess(rms,target_rms)

# @ddt
# class TestGravityB1Final(unittest.TestCase):
#     @classmethod
#     def setUpClass(cls):
#         cls.csm = CSM(dOmega0=0.32,dLambda=0.68,dSigma8=0.83,ns=0.96)
#         cls.msr = MSR()
#         cls.time = cls.msr.Load('b1-final-np-asym-k1.std')
#         cls.a = np.load('b1-final-np-asym-k1.acc.npy')
#         cls.maga = np.linalg.norm(cls.a,axis=1)
#         cls.msr.setParameters(bPeriodic=False,bEwald=False,nReplicas=None)

#     @data([0.70,0.001],)
#     @unpack
#     def testGravityBasic(self,theta,target_rms):
#         #print(self.msr.parm)
#         self.msr.DomainDecomp()
#         self.msr.BuildTree()
#         self.msr.Gravity(time=self.time,theta=theta) 
#         self.msr.Reorder()
#         a=self.msr.GetArray(field=2,time=self.time)
#         relerr = np.linalg.norm(a-self.a,axis=1) / self.maga
#         rms = np.std(relerr)
#         print('rms',rms)
#         self.assertLess(rms,target_rms)

if __name__ == '__main__':
    print('Running test')
    unittest.main(verbosity=2,testRunner=xmlrunner.XMLTestRunner(output='test-reports'))
