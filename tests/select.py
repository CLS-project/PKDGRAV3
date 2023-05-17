from __future__ import division
import os
from math import sqrt,pi
import numpy as np
import unittest
from ddt import ddt, data, file_data, unpack
import xmlrunner
from MASTER import MSR
from CSM import CSM

@ddt
@unittest.skipIf(not os.path.isfile('b0-final.std'), "missing b0-final.std")
#@unittest.skipIf(not os.path.isfile('b0-final-np-asym-k1.acc.npy'), "missing b0-final-np-asym-k1.acc.npy")
class TestSelect(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.csm = CSM(dOmega0=0.32,dLambda=0.68,dSigma8=0.83,ns=0.96)
        cls.msr = MSR()
        cls.time = cls.msr.Load('b0-final.std')

    # Spheres of different position and radius
    @data([[0.1,0.1,0.1],0.05,52],[[0.1,0.1,0.1],0.1,1700],)
    @unpack
    def testSelectSphere(self,center,radius,count):
        n = self.msr.MarkSphere(center=center,radius=radius)
        self.assertEqual(n,count)

    # Boxes
    @data([[0.1,0.1,0.1],[0.05,0.05,0.05],113],[[0.1,0.1,0.1],[0.1,0.1,0.1],17014],)
    @unpack
    def testSelectBox(self,center,size,count):
        n = self.msr.MarkBox(center=center,size=size)
        self.assertEqual(n,count)

    # Box and Sphere (overlapping and not)
    @data([[0.1,0.1,0.1],[0.05,0.05,0.05],[0.1,0.1,0.1],0.1,1700],[[0.1,0.1,0.1],[0.05,0.05,0.05],[-0.1,-0.1,-0.1],0.05,2339],)
    @unpack
    def testSelectBoxSphere(self,bcenter,bsize,scenter,sradius,count):
        n = self.msr.MarkBox(center=bcenter,size=bsize)
        n = self.msr.MarkSphere(center=scenter,radius=sradius,clearIfFalse=False)
        self.assertEqual(n,count)

    # Cylinder
    @data([[-0.1,0.1,0.1],[0.1,0.1,0.1],0.1,8195],[[-0.1,-0.1,-0.1],[0.1,0.1,0.05],0.15,142807],)
    @unpack
    def testSelectCylinder(self,point1,point2,radius,count):
        n = self.msr.MarkCylinder(point1=point1,point2=point2,radius=radius)
        self.assertEqual(n,count)



if __name__ == '__main__':
    print('Running test')
    unittest.main(verbosity=2,testRunner=xmlrunner.XMLTestRunner(output='test-reports'))
