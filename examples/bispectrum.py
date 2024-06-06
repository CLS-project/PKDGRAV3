# This is the prototype bispectrum analysis hook for pkdgrav3
import numpy as np
import os, pickle
from PKDGRAV import add_analysis
from accuracy import classic_theta_switch,classic_replicas_switch

# Class to perform analysis. It must be "callable".
class Bispectrum:

    from math import pi
    grid = 0
    order = 3
    def __init__(self,grid,nRadialBins,nThetaBins,achOutName,dBoxSize):
        self.grid = grid
        self.nRadialBins = nRadialBins
        self.nThetaBins  = nThetaBins
        self.achOutName = achOutName
        self.dBoxSize = dBoxSize

    def __call__(self,msr,step,time,**kwargs):
        self.powerspec(msr,step,time,**kwargs)
        self.bispec(msr,step,time,**kwargs)

    def ephemeral(self,msr,**kwargs):
        return msr.grid_ephemeral(self.grid,count=3+self.nThetaBins)

    def powerspec(self,msr,step,time,**kwargs):
        out_dir = 'powerspec'
        self.create_folder(out_dir) 

        msr.grid_create(self.grid)
        msr.assign_mass(grid_index=0)
        msr.density_contrast(grid_index=0)
        if True: # Interlace
            msr.assign_mass(grid_index=1,delta=0.5)
            msr.density_contrast(grid_index=1)
            msr.grid_interlace(target_grid_index=0,source_grid_index=1)
        msr.window_correction(grid_index=0)
        K,PK,N,*rest = msr.grid_bin_k(grid_index=0,bins=self.grid//2)
        msr.grid_delete()

        L = self.dBoxSize
        k_factor = 2.0 * self.pi / L
        pk_factor = L**3
        name = out_dir+'/{name}.{step:05d}.pk'.format(name=self.achOutName,step=step)
        with open(name,'w') as f:
            for k,pk,n in zip(K,PK,N):
                if n>0: f.write('{} {} {}\n'.format(k*k_factor,pk*pk_factor,n))

    def bispec(self,msr,step,time,**kwargs):
        ngrid = self.grid
        out_dir = 'bispec'
        self.create_folder(out_dir) 

        msr.grid_create(self.grid)
        msr.assign_mass(grid_index=0,order=self.order)
        # msr.grid_write(grid_index=0,filename="test.grid")
        msr.density_contrast(grid_index=0)
        if True: # Interlace
            msr.assign_mass(grid_index=1,delta=0.5,order=self.order)
            msr.density_contrast(grid_index=1)
            msr.grid_interlace(target_grid_index=0,source_grid_index=1)
        msr.window_correction(grid_index=0,order=self.order)

        L = self.dBoxSize
        print(f'BoxSize = {L} Mpc/h')
        k_factor = 2.0 * self.pi / L

        nRadialBins = self.nRadialBins
        kbin_edges = np.linspace(0,ngrid/4,int(nRadialBins)+1).astype(int) #[1.0,3.0,5.0,7.0]
        kbin_centr = np.sort([(kbin_edges[i]+kbin_edges[i+1])/2 for i in range(len(kbin_edges)-1)])
        
        nThetaBins = self.nThetaBins
        # theta = np.linspace(0, np.pi, nThetaBins)
        theta_edges = np.linspace(0, np.pi, nThetaBins+1)

        count = 0
        for p,k1 in enumerate(kbin_centr):
            for q,k2 in enumerate(kbin_centr[p:]):
                count += 1
                print(f'{count}/{int(len(kbin_centr)*(len(kbin_centr)+1)/2)}')
                # k3_centr = self.cosine_rule(k1, k2, np.pi-theta)
                # dk3 = np.gradient(kbin_centr)[0]
                # k3_dict = {i: (int(k-dk3/2),int(k+dk3/2)) for i,k in enumerate(k3_centr)}
                k3_edges = np.round(self.cosine_rule(k1, k2, np.pi-theta_edges))
                k3_centr = k3_edges[1:]/2+k3_edges[:-1]/2
                k3_dict = {i: (k3_edges[i],k3_edges[i+1]) for i,k in enumerate(k3_centr)}

                msr.bispectrum_normalize(kbin_edges[p],kbin_edges[p+1],target_grid_index=1)
                msr.bispectrum_normalize(kbin_edges[q],kbin_edges[q+1],target_grid_index=2)
                msr.bispectrum_select(   kbin_edges[p],kbin_edges[p+1],target_grid_index=1,source_grid_index=0)
                msr.bispectrum_select(   kbin_edges[q],kbin_edges[q+1],target_grid_index=2,source_grid_index=0)
                msr.bispectrum_normalize(kbin_edges[p],kbin_edges[p+1],target_grid_index=4)
                msr.bispectrum_normalize(kbin_edges[q],kbin_edges[q+1],target_grid_index=5)
                
                Bk12 = []
                out_file = out_dir+'/{0:}.{1:05d}.nGrid{2:}.k1_{3:.2f}_k2_{4:.2f}.bispec'.format(self.achOutName,step,ngrid,k1*k_factor,k2*k_factor)
                for r,k3 in enumerate(k3_centr):
                    # print(f'k1, k2, k3 = {k1*k_factor:.3f}, {k2*k_factor:.3f}, {k3*k_factor:.3f} h/Mpc')
                    print(f'k1, k2, k3 = {k1:}, {k2}, {k3:}')
                    kmin=k3_dict[r][0]
                    kmax=k3_dict[r][1]
                    if kmax-kmin<=1:
                        kmin -= 1
                        kmax += 1
                    if kmin<0: kmin = 0 
                    if kmax>ngrid: kmax = ngrid
                    msr.bispectrum_normalize(kmin,kmax,target_grid_index=3)
                    normalization = msr.bispectrum_calculate(4,5,3)
                    msr.bispectrum_select(kmin,kmax,target_grid_index=3,source_grid_index=0)
                    result = msr.bispectrum_calculate(1,2,3)
                    print('result, normalization =', result, normalization)
                    out = result/normalization
                    print('result/normalization =', out)
                    out = out*L**6
                    print('result/normalization*(V^2) =', out)
                    Bk12.append([k1*k_factor,k2*k_factor,k2*k_factor,out])
                
                np.savetxt(out_file, Bk12)

        msr.grid_delete()
        print('...done')

    def cosine_rule(self, k1, k2, alpha):
        k3 = np.sqrt(k1**2 + k2**2 + 2*k1*k2*np.cos(alpha))
        return k3
    
    def create_folder(self, grid_dir):
        if not os.path.exists(grid_dir):
            os.makedirs(grid_dir)
            print(f"Directory '{grid_dir}' created.")

# Here we add our analysis hook by contructing an instance of the object. We could have also just
# passed a function, but an object instance can be useful if we need state information like grid.
# We also need to specify how much memory we need by passing an object created with ephemeral().
# Without parameters this object indicates that we need no additional memory.
# Here we are saying that we need 50 grids of nGrid^3 plus one extra for the source.
_nGridBispectrum = 512 #64
_nRadialBins     = 16  #10
_nThetaBins      = 16

achOutName      = "example"

# Initial Condition
dBoxSize        = 256       # Mpc/h
nGrid           = 128        # Simulation has nGrid^3 particles
iLPT            = 2             # LPT order for IC
iSeed           = 314159265     # Random seed
dRedFrom        = 150         # Starting redshift

# Cosmology
achTfFile       = "euclid_z0_transfer_combined.dat"
h               = 0.67
dOmega0         = 0.32
dLambda         = 0.68
dSigma8         = 0.83
dSpectral       = 0.96

iStartStep      = 0
nSteps          = 100
dRedTo          = 0.0

# Cosmological Simulation
bComove         = True          # Use comoving coordinates
bPeriodic       = True          # with a periodic box
bEwald          = True          # enable Ewald periodic boundaries

# Logging/Output
iOutInterval    = 10
#iCheckInterval = 10
bDoDensity      = False
bVDetails       = True

bOverwrite      = True
bParaRead       = True          # Read in parallel
bParaWrite      = False         # Write in parallel (does not work on all file systems)
#nParaRead      = 8             # Limit number of simultaneous readers to this
#nParaWrite     = 8             # Limit number of simultaneous writers to this

# Accuracy Parameters
bEpsAccStep     = True          # Choose eps/a timestep criteria
dTheta          = classic_theta_switch()        # 0.40, 0.55, 0.70 switch
nReplicas       = classic_replicas_switch()     # 1 if theta > 0.52 otherwise 2

# Memory and performance
bMemUnordered   = True          # iOrder replaced by potential and group id
bNewKDK         = True          # No accelerations in the particle, dual tree possible

add_analysis(Bispectrum(grid=_nGridBispectrum,nRadialBins=_nRadialBins,nThetaBins=_nThetaBins,achOutName=achOutName,dBoxSize=dBoxSize))
