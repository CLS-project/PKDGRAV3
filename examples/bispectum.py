# This is the prototype bispectrum analysis hook for pkdgrav3

from MASTER import MSR # Once MSR is imported, simulation mode is no longer automatically entered
from math import pi

msr = MSR()

# Read in the parameters from another file. We will run the standard "cosmology.par" simulation
from importlib.machinery import SourceFileLoader
params = vars(SourceFileLoader('cosmology', 'pk.par').load_module())
msr.setParameters(**params)

# Class to perform analysis. It must be "callable".
class bispectrum:
    grid = 0
    order = 4
    def __init__(self,grid):
        self.grid = grid
    def __call__(self,msr,step,time,a,**kwargs):
        print('Calculating bispectrum')
        msr.grid_prepare(self.grid)
        msr.assign_mass(target=0,order=self.order)
        msr.density_contrast(target=0)
        if True: # Interlace
            msr.assign_mass(target=1,delta=0.5,order=self.order)
            msr.density_contrast(target=1)
            msr.interlace(target=0,source=1)
        msr.window_correction(target=0,order=self.order)
        msr.bispectrum_select(target=1,source=0,kmin=1.0,kmax=3.0)
        msr.bispectrum_select(target=2,source=0,kmin=3.0,kmax=5.0)
        msr.bispectrum_select(target=3,source=0,kmin=5.0,kmax=7.0)
        result = msr.bispectrum_calculate(1,2,3)
        print(result)
        msr.grid_release()

# Here we add our analysis hook by contructing an instance of the object. We could have also just
# passed a function, but an object instance can be useful if we need state information like grid.
# We also need to specify how much memory we need by passing an object created with ephemeral().
# Without parameters this object indicates that we need no additional memory.
# Here we are saying that we need 50 grids of nGrid^3 plus one extra for the source.
nGridBispectrum=64
nRadialBins=50
msr.add_analysis(bispectrum(grid=nGridBispectrum),msr.ephemeral(grid=nGridBispectrum,count=1 + nRadialBins))

# Since simulation mode is not automatically invoked, we need to do it manually.
msr.simulate()
