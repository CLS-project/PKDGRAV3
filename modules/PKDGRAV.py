# distutils: language = c++
# cython: always_allow_keywords=True

import cython
# from cython.cimports.cpython import array
# import array
import numpy as np
from cosmology import Cosmology

def set_parameters(**kwargs):
    msr0.parameters.update(kwargs,True)

def restart(arguments,specified,species,classes,n,name,step,steps,time,delta,E,U,Utime):
    ndark = cython.declare(cython.size_t,species[FIO_SPECIES.FIO_SPECIES_DARK])
    nsph  = cython.declare(cython.size_t,species[FIO_SPECIES.FIO_SPECIES_SPH])
    nstar = cython.declare(cython.size_t,species[FIO_SPECIES.FIO_SPECIES_STAR])
    nbh   = cython.declare(cython.size_t,species[FIO_SPECIES.FIO_SPECIES_BH])
    aClasses = new_partclass_vector()
    for r in classes:
        spec = cython.declare(cython.int,  r[0])
        mass = cython.declare(cython.float,r[1])
        soft = cython.declare(cython.float,r[2])
        imat = cython.declare(cython.int,  r[3])
        aClasses.push_back(PARTCLASS(FIO_SPECIES(spec),mass,soft,imat))
    msr0.Restart(n,name.encode('UTF-8'),
        step,steps,time,delta,ndark,nsph,nstar,nbh,
        E,U,Utime,aClasses,arguments,specified)

def generate_ic(cosmology : Cosmology,*,grid : int,seed : int,z : float,L : float,
                order : int = None, fixed_amplitude : bool = False, phase_pi : float = 0, **kwargs) -> float:
    """
    Generate initial conditions for a cosmological simulation.
    
    :param Cosmology cosmology: cosmology
    :param integer grid: grid size of the initial conditions
    :param integer seed: random seed
    :param number z: starting redshift
    :param number L: length unit of the box
    :param integer order: IC order, 1=Zeldovich, 2=2LPT
    :param Boolean fixed_amplitude: use fixed amplitude for the power spectrum
    :param number phase_pi: phase of the initial conditions (in units of :math:`\\pi` radians, normally 0 or 1)
    :return: time
    """
    msr0.parameters.set(msr0.parameters.str_bFixedAmpIC,fixed_amplitude)
    msr0.parameters.set(msr0.parameters.str_dFixedAmpPhasePI,phase_pi)
    if order is None: pass
    elif order == 1:  msr0.parameters.set(msr0.parameters.str_b2LPT,False)
    elif order == 2:  msr0.parameters.set(msr0.parameters.str_b2LPT,True)
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

def domain_decompose(rung=0):
    """
    Particles are ordered spatially across all nodes and cores using
    the Orthagonal Recursive Bisection (ORB) method.

    :param integer rung: balance work based on this rung
    """
    msr0.DomainDecomp(rung)

def build_tree(ewald=False):
    """
    Builds a tree in each domain

    :param Boolean ewald: construct a moment for Ewald if true
    """
    msr0.BuildTree(ewald)

def reorder():
    """
    Returns particles to their original order (usually done for output).
    If particles are unordered then this function has no effect.
    """
    msr0.Reorder()


def gravity(time=0.0,delta=0.0,theta=0.7,rung=0,ewald=None,step=0.0,kick_close=True,kick_open=True, only_marked=False):
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
    msr0.Hostname()
    dTime = msr0.LoadOrGenerateIC()
    if dTime >= 0:
        msr0.Simulate(dTime)

def measure_pk(grid,bins=0,a=1.0,interlace=True,order=4,L=1.0):
    """
    Measure the Power spectrum P(k) for the box.

    :param integer grid: grid size for the mass assignment
    :param integer bins: number of bins for P(k), defaults to half the grid size
    :param number a: expansion factor
    :param Boolean interlace: use interlacing to reduce grid aliasing
    :param integer order: mass assignment order, 1=NGP, 2=CIC, 3=TSC, 4=PCS
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

def fof(tau,minmembers=10):
    """
    Friends of friends (fof) group finding

    :param number tau: linking length
    :param integer minmembers: minimum group size (in particles)
    """
    msr0.NewFof(tau,minmembers)
    msr0.GroupStats()

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
