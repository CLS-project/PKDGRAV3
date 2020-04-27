# This is an example of an analysis hook. It does essentially the same as the built in P(k) measurement.

from MASTER import MSR # Once MSR is imported, simulation mode is no longer automatically entered
from math import pi

msr = MSR()

# Read in the parameters from another file. We will run the standard "cosmology.par" simulation
from importlib.machinery import SourceFileLoader
params = vars(SourceFileLoader('cosmology', 'cosmology.par').load_module())
msr.setParameters(**params)

# Here is another way
# Normal pkdgrav3 parameters should be set in the function.
# This way we avoid naming conflicts.
#def getParameters():
#    bEpsAccStep = True
#    ... and all others
#    return(locals())
#msr.setParameters(**getParameters())

# Class to perform analysis. It must be "callable".
class measurepk:
    grid = 0
    def __init__(self,grid):
        self.grid = grid
        print('adding P(k) analysis with grid {grid}'.format(grid=grid))
    def __call__(self,msr,step,time,a,**kwargs):
        msr.grid_prepare(self.grid)
        msr.assign_mass(target=0)
        msr.density_contrast(target=0)
        if True: # Interlace
            msr.assign_mass(target=1,delta=0.5)
            msr.density_contrast(target=1)
            msr.interlace(target=0,source=1)
        msr.window_correction(target=0)
        K,PK,N,*rest = msr.grid_bin_k(source=0,bins=self.grid//2)
        msr.grid_release()

        L = msr.parm.dBoxSize
        k_factor = 2.0 * pi / L
        pk_factor = L**3
        name='{name}.{step:05d}.pk'.format(name=msr.parm.achOutName,step=step)
        with open(name,'w') as f:
            for k,pk,n in zip(K,PK,N):
                if n>0: f.write('{} {} {}\n'.format(k*k_factor,pk*pk_factor,n))

# Here we add our analysis hook by contructing an instance of the object. We could have also just
# passed a function, but an object instance can be useful if we need state information like grid.
# We also need to specify how much memory we need by passing an object created with ephemeral().
# Without parameters this object indicates that we need no additional memory.
# Here we are saying that we need two grids of nGrid^3.
msr.add_analysis(measurepk(grid=params['nGrid']),msr.ephemeral(grid=params['nGrid'],count=2))

# Since simulation mode is not automatically invoked, we need to do it manually.
msr.simulate()
