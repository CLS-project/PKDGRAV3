# distutils: language = c++
# cython: always_allow_keywords=True

import cython
import sys
# from cython.cimports.cpython import array
# import array
import numpy as np
from cosmology import Cosmology

if 'sphinx' not in sys.modules:
    import ephemeral

def arguments():
    return msr0.parameters.arguments()

def set_parameters(**kwargs):
    if not msr0.parameters.update(kwargs,False):
        raise ValueError("invalid parameter")

def add_analysis(callable,memory=None):
    """
    Add an analysis function to the simulation.

    :param callable: the analysis function
    :param memory: the ephemeral memory required by the analysis function

    If memory is not specified then an attempt is made to use the ephemeral
    method of the callable to get the memory requirements if it exists.
    Otherwise it is assumed that the callable does not require ephemeral memory.
    """
    if memory is None:
        if hasattr(callable, 'ephemeral'):
            memory = callable.ephemeral(__import__('PKDGRAV'))
        else:
            memory = ephemeral.PyEphemeralMemory(0,0)
    if not isinstance(memory,ephemeral.PyEphemeralMemory):
        raise ValueError("invalid ephemeral memory")
    msr0.addAnalysis(callable,memory.per_particle,memory.per_process)

def restore(filename,species=None,classes=None,step=None,steps=None,time=None,delta=None,E=None,U=None,Utime=None,**kwargs):
    """
    Restore a simulation from a file.

    :param str filename: the name of the file
    :param species: species counts (optional)
    :param classes: particle classes (optional)
    :param integer step: current step (optional)
    :param integer steps: total steps (optional)
    :param number time: simulation time (optional)
    :param number delta: time step delta (optional)
    :param number E: total energy (optional)
    :param number U: potential energy (optional)
    :param number Utime: time of potential energy calculation (optional)

    The species and classes parameters are lists of tuples. Each tuple
    contains the species number, mass, softening, and imat values.

    The kwargs are additional parameters to set.
    """
    msr0.Restart(filename.encode('UTF-8'),kwargs,species,classes,step,steps,time,delta,E,U,Utime)

def generate_ic(cosmology : Cosmology,*,grid : int,seed : int,z : float,L : float,
                order : int = None, fixed_amplitude : bool = False, phase_pi : float = 0, **kwargs) -> float:
    """
    Generate initial conditions for a cosmological simulation.
    
    :param Cosmology cosmology: cosmology
    :param integer grid: grid size of the initial conditions
    :param integer seed: random seed
    :param number z: starting redshift
    :param number L: length unit of the box
    :param integer order: IC order, 1=Zeldovich, 2=2LPT, 3=3LPT
    :param Boolean fixed_amplitude: use fixed amplitude for the power spectrum
    :param number phase_pi: phase of the initial conditions (in units of :math:`\\pi` radians, normally 0 or 1)
    :return: time
    """
    msr0.parameters.set(msr0.parameters.str_bFixedAmpIC,fixed_amplitude)
    msr0.parameters.set(msr0.parameters.str_dFixedAmpPhasePI,phase_pi)
    if order is None: pass
    elif order == 1: msr0.parameters.set(msr0.parameters.str_iLPT, 1)
    elif order == 2: msr0.parameters.set(msr0.parameters.str_iLPT, 2)
    elif order == 3: msr0.parameters.set(msr0.parameters.str_iLPT, 3)
    else:             raise ValueError("invalid IC order")
    set_parameters(**kwargs)
    return msr0.GenerateIC(grid,seed,z,L,cosmology._csm)

def load(filename,**kwargs):
    """
    Read particles from an input file.

    :param str filename: the name of the file
    :return: time
    :rtype: number
    """
    set_parameters(**kwargs)
    return msr0.Read(filename.encode('UTF-8'))

def save(filename,time=1.0):
    """
    Save particles to a file.

    :param str filename: the name of the file
    :param number time: simulation time
    """
    return msr0.Write(filename.encode('UTF-8'),time,False)

def load_checkpoint(filename,species=None,classes=None,step=None,steps=None,time=None,delta=None,E=None,U=None,Utime=None,**kwargs):
    """
    Read particles from a checkpoint file.

    :param str filename: the name of the file
    :return: time
    :rtype: number
    """
    return read_checkpoint(filename,kwargs,species,classes,step,steps,time,delta,E,U,Utime)

def domain_decompose(rung=0):
    """
    Particles are ordered spatially across all nodes and cores using
    the Orthagonal Recursive Bisection (ORB) method.

    :param integer rung: balance work based on this rung
    """
    msr0.DomainDecomp(rung)

def build_tree(ewald=None):
    """
    Builds a tree in each domain

    :param Boolean ewald: construct a moment for Ewald if true
    """
    bEwald = msr0.parameters.get_bEwald() if ewald is None else ewald
    msr0.BuildTree(bEwald)

def reorder():
    """
    Returns particles to their original order (usually done for output).
    If particles are unordered then this function has no effect.
    """
    msr0.Reorder()


def gravity(time=0.0,delta=0.0,theta=0.7,rung=0,ewald=None,step=0.0,kick_close=True,kick_open=False, only_marked=False):
    bEwald = msr0.parameters.get_bEwald() if ewald is None else ewald
    bGravStep = msr0.parameters.get_bGravStep()
    nPartRhoLoc = msr0.parameters.get_nPartRhoLoc()
    iTimeStepCrit = msr0.parameters.get_iTimeStepCrit()
    r= msr0.Gravity(rung,63,1,
                        3 if only_marked else 0,
                        time,delta,step,theta,kick_close,kick_open,bEwald,
                        bGravStep,nPartRhoLoc,iTimeStepCrit)
    # return r

def simulate(**kwargs):
    """
    Directly enter simulation mode. Normally simulation mode is disabled
    if PKDGRAV has been imported. Calling this function proceeds normally.
    """
    set_parameters(**kwargs)
    msr0.ValidateParameters()
    msr0.Hostname()
    dTime = msr0.LoadOrGenerateIC()
    if dTime >= 0:
        msr0.Simulate(dTime)

def measure_pk(grid,bins=0,a=1.0,interlace=True,order=3,L=1.0):
    """
    Measure the Power spectrum P(k) for the box.

    :param integer grid: grid size for the mass assignment
    :param integer bins: number of bins for P(k), defaults to half the grid size
    :param number a: expansion factor
    :param Boolean interlace: use interlacing to reduce grid aliasing
    :param integer order: mass assignment order, 0=NGP, 1=CIC, 2=TSC, 3=PCS
    :param number L: length unit of the box to convert k and P(k) to physical units
    :return: k, P(k), N(k), Pall(k)
    :rtype: tuple of numpy arrays
    """
    from math import pi
    if bins==0: bins=grid//2
    (npk,k,pk,lpk) = MeasurePk(order,interlace,grid,a,bins)
    k *= 2.0 * pi / L
    pk *= L**3
    lpk *= L**3
    return (k,pk,npk,lpk)

def grid_bin_k(bins,grid_index):
    """
    Bin the k-space grid

    :param integer bins: number of bins
    :param integer grid_index: which grid number to use
    """
    return GridBinK(bins,grid_index)

def grid_ephemeral(grid,count=1):
    """
    Return the Ephemeral memory required for the grid

    :param integer grid: grid size for the mass assignment
    :param integer count: number of grids
    """
    result = msr0.EphemeralMemoryGrid(grid,count)
    return ephemeral.PyEphemeralMemory(result.per_particle,result.per_process)

def grid_create(grid):
    """
    Create a grid for the mass assignment

    :param integer grid: grid size for the mass assignment
    """
    msr0.GridCreateFFT(grid)

def grid_delete():
    """
    Delete the grid for the mass assignment
    """
    msr0.GridDeleteFFT()

def grid_write(filename,k=False,grid_index=0):
    """
    Write the grid to a file

    :param str filename: the name of the file
    :param Boolean k: write the k-space grid
    :param integer grid_index: grid index
    :param Boolean parallel: number of parallel tasks
    """
    msr0.OutputGrid(filename.encode('UTF-8'), k, grid_index, 1)

def assign_mass(order=3,grid_index=0,delta=0.0,fold=1):
    """
    Assign mass to the grid

    :param integer order: mass assignment order, 0=NGP, 1=CIC, 2=TSC, 3=PCS
    :param integer grid_index: which grid number to use
    :param number delta: grid shift (normally 0.0 or 0.5)
    :param number fold: number of times to fold
    """
    msr0.AssignMass(order,grid_index,delta,fold)

def density_contrast(grid_index=0,k=True):
    """
    Compute the density contrast

    :param integer grid_index: which grid number to use
    """
    msr0.DensityContrast(grid_index,k)

def grid_interlace(target_grid_index=0,source_grid_index=0):
    """
    Interlace the grid

    :param integer target_grid_index: which grid number to use
    :param integer source_grid_index: which grid number to use
    """
    msr0.Interlace(target_grid_index,source_grid_index)

def window_correction(grid_index=0,order=3):
    """
    Apply window correction to the grid

    :param integer grid_index: which grid number to use
    :param integer order: mass assignment order, 0=NGP, 1=CIC, 2=TSC, 3=PCS
    """
    msr0.WindowCorrection(order,grid_index) 

def add_linear_signal(seed,L,a,fixed=False,phase=0.0,grid_index=0):
    """
    Add a linear signal to the grid

    :param integer seed: random seed
    :param number L: length unit of the box
    :param number a: expansion factor
    :param Boolean fixed: use fixed amplitude for the power spectrum
    :param number phase: phase of the initial conditions (in units of :math:`\\pi` radians, normally 0 or 1)
    :param integer grid_index: which grid number to use (default 0)
    """
    msr0.AddLinearSignal(grid_index,seed,L,a,fixed,phase)


def bispectrum_select(k_min,k_max,target_grid_index=0,source_grid_index=1):
    """
    Select the bispectrum

    :param integer target_grid_index: which grid number to use
    :param integer source_grid_index: which grid number to use
    :param number k_min: minimum k value
    :param number k_max: maximum k value
    """
    msr0.BispectrumSelect(target_grid_index,source_grid_index,k_min,k_max)

def bispectrum_normalize(k_min,k_max,target_grid_index=0):
    """
    Normalize the bispectrum

    :param integer target_grid_index: which grid number to use
    :param number k_min: minimum k value
    :param number k_max: maximum k value
    """
    msr0.BispectrumSelect(target_grid_index,-1,k_min,k_max)

def bispectrum_calculate(grid_index0,grid_index1,grid_index2):
    """
    Calculate the bispectrum

    :param integer grid_index0: which grid number to use
    :param integer grid_index1: which grid number to use
    :param integer grid_index2: which grid number to use
    """
    return msr0.BispectrumCalculate(grid_index0,grid_index1,grid_index2)

def fof(tau,minmembers=10):
    """
    Friends of friends (fof) group finding

    :param number tau: linking length
    :param integer minmembers: minimum group size (in particles)
    """
    msr0.NewFof(tau,minmembers)
    msr0.GroupStats()

def smooth(type,n=32,time=1.0,delta=0.0,symmetric=False,resmooth=False):
    """
    Smooths the density field with a given kernel

    :param integer type: smoothing kernel type
    :param integer n: smoothing kernel size
    :param number time: simulation time
    :param number delta: time step delta
    :param Boolean symmetric: use symmetric smoothing

    Values for the smoothing kernel type are:

    * SMOOTH_TYPE_DENSITY
    * SMOOTH_TYPE_F1
    * SMOOTH_TYPE_M3
    * SMOOTH_TYPE_GRADIENT_M3
    * SMOOTH_TYPE_HOP_LINK
    * SMOOTH_TYPE_BALL
    * SMOOTH_TYPE_PRINTNN
    * SMOOTH_TYPE_HYDRO_DENSITY
    * SMOOTH_TYPE_HYDRO_DENSITY_FINAL
    * SMOOTH_TYPE_HYDRO_GRADIENT
    * SMOOTH_TYPE_HYDRO_FLUX
    * SMOOTH_TYPE_HYDRO_STEP
    * SMOOTH_TYPE_HYDRO_FLUX_VEC
    * SMOOTH_TYPE_SN_FEEDBACK
    * SMOOTH_TYPE_BH_MERGER
    * SMOOTH_TYPE_BH_DRIFT
    * SMOOTH_TYPE_BH_STEP
    * SMOOTH_TYPE_CHEM_ENRICHMENT
    """
    if resmooth:
        msr0.ReSmooth(time,delta,type,symmetric)
    else:
        msr0.Smooth(time,delta,type,symmetric,n)

def get_array(field,time=1.0,marked=False):
    """
    Retrieves an array with requested field.

    :param number field: the field to retrieve. Values are:

    * FIELD_POSITION
    * FIELD_ACCELERATION
    * FIELD_VELOCITY
    * FIELD_POTENTIAL
    * FIELD_GROUP
    * FIELD_MASS
    * FIELD_SOFTENING
    * FIELD_DENSITY
    * FIELD_BALL
    * FIELD_PARTICLE_ID
    * FIELD_GLOBAL_GID

    :param number time: simulation time
    :param Boolean marked: retrieve only marked particles
    """
    N = np.array([msr0.N,1],dtype=np.uint64)
    T = np.float32
    if marked: N[0] = msr0.CountSelected()
    if field == FIELD_POSITION:
        N[1] = 3
        T = np.float64
    elif field == FIELD_ACCELERATION:
        N[1] = 3
    elif field == FIELD_VELOCITY:
        N[1] = 3
    elif field == FIELD_POTENTIAL:
        pass
    elif field == FIELD_GROUP:
        T = np.int32
    elif field == FIELD_MASS:
        pass
    elif field == FIELD_SOFTENING:
        pass
    elif field == FIELD_DENSITY:
        pass
    elif field == FIELD_BALL:
        pass
    elif field == FIELD_PARTICLE_ID:
        T = np.uint64
    elif field == FIELD_GLOBAL_GID:
        T = np.uint64
    else:
        raise ValueError("invalid array requested")
    a = np.zeros(N,dtype=T)
    if   T == np.float32: v = a2f2(a)
    elif T == np.float64: v = a2d2(a)
    else:                 v = a2u2(a)
    msr0.RecvArray(v,field,N[1]*a.itemsize,time,marked)
    if N[1] == 1: a = np.reshape(a,(N[0]))
    return a

def write_array(filename,field):
    """
    Writes an array to a file.

    :param str filename: the name of the file
    :param number field: the field to write. Values are:

    * OUT_DENSITY_ARRAY
    * OUT_POT_ARRAY
    * OUT_AMAG_ARRAY
    * OUT_IMASS_ARRAY
    * OUT_RUNG_ARRAY
    * OUT_DIVV_ARRAY
    * OUT_VELDISP2_ARRAY
    * OUT_VELDISP_ARRAY
    * OUT_PHASEDENS_ARRAY
    * OUT_SOFT_ARRAY
    * OUT_POS_VECTOR
    * OUT_VEL_VECTOR
    * OUT_ACCEL_VECTOR
    * OUT_MEANVEL_VECTOR
    * OUT_IORDER_ARRAY
    * OUT_C_ARRAY
    * OUT_HSPH_ARRAY
    * OUT_RUNGDEST_ARRAY
    * OUT_MARKED_ARRAY
    * OUT_CACHEFLUX_ARRAY
    * OUT_CACHECOLL_ARRAY
    * OUT_AVOIDEDFLUXES_ARRAY
    * OUT_COMPUTEDFLUXES_ARRAY
    * OUT_HOP_STATS
    * OUT_GROUP_ARRAY
    * OUT_GLOBALGID_ARRAY
    * OUT_BALL_ARRAY
    * OUT_PSGROUP_ARRAY
    * OUT_PSGROUP_STATS
    """
    msr0.OutASCII(filename.encode('UTF-8'),field,3 if field in [OUT_POS_VECTOR,OUT_VEL_VECTOR,OUT_MEANVEL_VECTOR,OUT_ACCEL_VECTOR] else 1,0)

def mark_box(center,apothem,set_if_true=1,clear_if_false=1):
    """
    Mark particles inside a given box

    :param number[3] center: center coordinate
    :param number[3] apothem: distance from center to edge
    :param integer set_if_true: mark the particle if it is inside
    :param integer clear_if_flase: unmark the particle if it is outside
    :return: number of particles marked
    :rtype: integer
    """
    return msr0.SelBox(TinyVector[double,BLITZ3](center[0],center[1],center[2]),
                TinyVector[double,BLITZ3](apothem[0],  apothem[1],  apothem[2]),
                set_if_true,clear_if_false)

def mark_sphere(center,radius,set_if_true=1,clear_if_false=1):
    """
    Mark particles inside a given sphere

    :param number[3] center: center coordinate
    :param number radius: distance from center to edge
    :param integer set_if_true: mark the particle if it is inside
    :param integer clear_if_flase: unmark the particle if it is outside
    :return: number of particles marked
    :rtype: integer
    """
    return msr0.SelSphere(TinyVector[double,BLITZ3](center[0],center[1],center[2]),
                   radius,set_if_true,clear_if_false)

def mark_cylinder(point1,point2,radius,set_if_true=1,clear_if_false=1):
    """
    Mark particles inside a cylinder

    :param number[3] point1: center of the first end of the cylinder
    :param number[3] point2: center of the second end of the cylinder
    :param number radius: radius of the cylinder
    :param integer set_if_true: mark the particle if it is inside
    :param integer clear_if_flase: unmark the particle if it is outside
    :return: number of particles marked
    :rtype: integer
    """
    return msr0.SelCylinder(TinyVector[double,BLITZ3](point1[0],point1[1],point1[2]),
                     TinyVector[double,BLITZ3](point2[0],point2[1],point2[2]),
                     radius,set_if_true,clear_if_false)

def mark_species(*,species,species_mask,set_if_true=1,clear_if_false=1):
    """
    Mark particles of a given species

    :param integer species: species number
    :param integer species_mask: species mask (one or more species)
    :param integer set_if_true: mark the particle if it is inside
    :param integer clear_if_flase: unmark the particle if it is outside
    :return: number of particles marked
    :rtype: integer
    """
    if species_mask is None:
        species_mask = 1 << species
    return msr0.SelSpecies(species_mask,set_if_true,clear_if_false)

def rs_load_ids(filename,append=False):
    """
    Load particle IDs from a file

    :param str filename: the name of the file
    :param Boolean bAppend: append to the existing IDs
    """
    msr0.RsLoadIds(filename.encode('UTF-8'),append)

def rs_halo_load_ids(filename,append=False):
    """
    Load halo IDs from a file

    :param str filename: the name of the file
    :param Boolean bAppend: append to the existing IDs
    """
    msr0.RsHaloLoadIds(filename.encode('UTF-8'),append)

def rs_save_ids(filename):
    """
    Save particle IDs to a file

    :param str filename: the name of the file
    """
    msr0.RsSaveIds(filename.encode('UTF-8'))

def rs_reorder_ids():
    """
    Reorder the IDs
    """
    msr0.RsReorderIds()

def rs_extract(filename):
    """
    Extract particles to a file

    :param str filename_template: the name of the file
    """
    msr0.RsExtract(filename.encode('UTF-8'))