from __future__ import division
import os
from math import sqrt,pi
import numpy as np
import unittest
import xmlrunner
import PKDGRAV as msr
from PKDGRAV import HYDRO_MODEL
from importlib.machinery import SourceFileLoader

@unittest.skipIf(not os.path.isfile('collision.par'), "missing collision.par")
@unittest.skipIf(not os.path.isfile('collision.std'), "missing collision.std")
@unittest.skipIf(not os.path.isfile('MANEOStable_iron.in'), "missing MANEOStable_iron.in")
@unittest.skipIf(not os.path.isfile('MANEOStable_fosterite.in'), "missing MANEOStable_fosterite.in")
@unittest.skipIf(not os.path.isfile('initial_rho.npy'), "missing initial_rho.npy")
@unittest.skipIf(not os.path.isfile('initial_fBall.npy'), "missing initial_fBall.npy")
class TestNewSPHRead(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.params = vars(SourceFileLoader('collision', 'collision.par').load_module())
        cls.params['bGasDoExtensiveILPTest'] = 0
        cls.rho = np.load('initial_rho.npy')
        cls.fBall = np.load('initial_fBall.npy')

    def testNewSPHRead(self):
        msr.set_parameters(hydro_model=HYDRO_MODEL.SPH,**self.params)
        msr.load(self.params['achInFile'])
        msr.reorder()
        # Get scalar values
        rho = msr.get_array(field=msr.FIELD_DENSITY,time=0.0)
        fBall = msr.get_array(field=msr.FIELD_BALL,time=0.0)
        relDiff_rho = np.abs((self.rho - rho) / self.rho)
        relDiff_fBall = np.abs((self.fBall - fBall) / self.fBall)
        print("Maximum relative density deviation = {}".format(np.max(relDiff_rho)))
        print("rms density deviation = {}".format(np.std(relDiff_rho)))
        print("Maximum relative fBall deviation = {}".format(np.max(relDiff_fBall)))
        print("rms fBall deviation = {}".format(np.std(relDiff_fBall)))
        self.assertLess(np.max(relDiff_rho),3e-4)
        self.assertLess(np.std(relDiff_rho),4e-7)
        self.assertLess(np.max(relDiff_fBall),6e-5)
        self.assertLess(np.std(relDiff_fBall),4e-6)

@unittest.skipIf(not os.path.isfile('collision.par'), "missing collision.par")
@unittest.skipIf(not os.path.isfile('collision.std'), "missing collision.std")
@unittest.skipIf(not os.path.isfile('MANEOStable_iron.in'), "missing MANEOStable_iron.in")
@unittest.skipIf(not os.path.isfile('MANEOStable_fosterite.in'), "missing MANEOStable_fosterite.in")
@unittest.skipIf(not os.path.isfile('final_pos.npy'), "missing final_rho.npy")
@unittest.skipIf(not os.path.isfile('final_vel.npy'), "missing final_rho.npy")
@unittest.skipIf(not os.path.isfile('final_acc.npy'), "missing final_acc.npy")
@unittest.skipIf(not os.path.isfile('final_rho.npy'), "missing final_rho.npy")
@unittest.skipIf(not os.path.isfile('final_fBall.npy'), "missing final_fBall.npy")
class TestNewSPHStep(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.params = vars(SourceFileLoader('collision', 'collision.par').load_module())
        cls.params['bGasDoExtensiveILPTest'] = 0
        cls.pos = np.load('final_pos.npy')
        cls.vel = np.load('final_vel.npy')
        cls.acc = np.load('final_acc.npy')
        cls.rho = np.load('final_rho.npy')
        cls.fBall = np.load('final_fBall.npy')

    def testNewSPHStep(self):
        msr.set_parameters(hydro_model=HYDRO_MODEL.SPH,**self.params)
        msr.simulate()
        msr.reorder()
        # Get vector values
        pos = msr.get_array(field=msr.FIELD_POSITION,time=0.0)
        vel = msr.get_array(field=msr.FIELD_VELOCITY,time=0.0)
        acc = msr.get_array(field=msr.FIELD_ACCELERATION,time=0.0)
        relDiff_pos = np.abs(np.linalg.norm(self.pos - pos, axis=1) / np.linalg.norm(self.pos, axis=1))
        relDiff_vel = np.abs(np.linalg.norm(self.vel - vel, axis=1) / np.linalg.norm(self.vel, axis=1))
        relDiff_acc = np.abs(np.linalg.norm(self.acc - acc, axis=1) / np.linalg.norm(self.acc, axis=1))
        # Get scalar values
        rho = msr.get_array(field=msr.FIELD_DENSITY,time=0.0)
        fBall = msr.get_array(field=msr.FIELD_BALL,time=0.0)
        relDiff_rho = np.abs((self.rho - rho) / self.rho)
        relDiff_fBall = np.abs((self.fBall - fBall) / self.fBall)
        print("Maximum relative position deviation = {}".format(np.max(relDiff_pos)))
        print("rms position deviation = {}".format(np.std(relDiff_pos)))
        print("Maximum relative velocity deviation = {}".format(np.max(relDiff_vel)))
        print("rms velocity deviation = {}".format(np.std(relDiff_vel)))
        print("Maximum relative acceleration deviation = {}".format(np.max(relDiff_acc)))
        print("rms acceleration deviation = {}".format(np.std(relDiff_acc)))
        print("Maximum relative density deviation = {}".format(np.max(relDiff_rho)))
        print("rms density deviation = {}".format(np.std(relDiff_rho)))
        print("Maximum relative fBall deviation = {}".format(np.max(relDiff_fBall)))
        print("rms fBall deviation = {}".format(np.std(relDiff_fBall)))
        self.assertLess(np.max(relDiff_pos),4e-7)
        self.assertLess(np.std(relDiff_pos),2e-9)
        self.assertLess(np.max(relDiff_vel),1e-4)
        self.assertLess(np.std(relDiff_vel),4e-7)
        self.assertLess(np.max(relDiff_acc),4e-2)
        self.assertLess(np.std(relDiff_acc),2e-4)
        self.assertLess(np.max(relDiff_rho),4e-4)
        self.assertLess(np.std(relDiff_rho),8e-7)
        self.assertLess(np.max(relDiff_fBall),5e-5)
        self.assertLess(np.std(relDiff_fBall),3e-6)

@unittest.skipIf(not os.path.isfile('collision_small.par'), "missing collision_small.par")
@unittest.skipIf(not os.path.isfile('collision_small.std'), "missing collision_small.std")
@unittest.skipIf(not os.path.isfile('MANEOStable_iron.in'), "missing MANEOStable_iron.in")
@unittest.skipIf(not os.path.isfile('MANEOStable_fosterite.in'), "missing MANEOStable_fosterite.in")
class TestNewSPHExtensive(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.params = vars(SourceFileLoader('collision_small', 'collision_small.par').load_module())
        cls.params['bGasDoExtensiveILPTest'] = 1

    def testNewSPHExtensive(self):
        msr.set_parameters(hydro_model=HYDRO_MODEL.SPH,**self.params)
        msr.simulate()
        self.assertLess(0.1,1.0)

if __name__ == '__main__':
    print('Running test')
    unittest.main(verbosity=2,testRunner=xmlrunner.XMLTestRunner(output='test-reports'))
