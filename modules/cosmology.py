# distutils: language = c++
# cython: always_allow_keywords=True

import cython
from math import sqrt,pi

class Cosmology:
    def __cinit__(self):
        self.initialize()
    def __dealloc__(self):
        self.finish()
    def __enter__(self):
        return self
    def __exit__(self, type, value, traceback):
        pass

    @property
    def omega_matter(self):
        return self._csm.val.dOmega0

    @property
    def omega_baryon(self):
        return self._csm.val.dOmegab

    @property
    def omega_lambda(self):
        return self._csm.val.dLambda

    @property
    def omega_radiation(self):
        return self._csm.val.dOmegaRad

    @property
    def omega_dark_energy(self):
        return self._csm.val.dOmegaDE

    @property
    def w0(self):
        return self._csm.val.w0

    @property
    def wa(self):
        return self._csm.val.wa

    def radiation_matter_equivalence(self):
        return cosmo.csmRadMatEquivalence(self._csm)
        
    def expansion_to_hubble(self,a):
        return cosmo.csmExp2Hub(self._csm,a)

    def time_to_hubble(self,time):
        return cosmo.csmTime2Hub(self._csm,time)

    def expansion_to_time(self,a):
        return cosmo.csmExp2Time(self._csm,a)

    def time_to_expansion(self,time):
        return cosmo.csmTime2Exp(self._csm,time)

    def comoving_drift_factor(self,time,delta):
        return cosmo.csmComoveDriftFac(self._csm, time, delta);

    def comoving_kick_factor(self,time,delta):
        return cosmo.csmComoveKickFac(self._csm, time, delta);

    def comoving_growth(self,a):
        return self._comoving_growth(a)

class SimpleCosmology(Cosmology):
    def __init__(self,*,omega_matter,omega_lambda=None,omega_baryon=0.0,omega_radiation=0.0,omega_dark_energy=0.0,
                        w0=-1.0,wa=0.0,running=0.0,pivot=0.05,comoving=True,sigma8=0.0,As=0.0,ns=0.0,h=1.0):
        self.initialize_simple(omega_matter,omega_lambda if omega_lambda is not None else 1.0 - omega_matter - omega_radiation - omega_dark_energy,
                               omega_baryon,omega_radiation,omega_dark_energy,w0,wa,running,pivot,comoving,sigma8,As,ns,h)

class ClassCosmology(Cosmology):
    def __init__(self,file,L=1.0,As=0.0,ns=0.0,linear_species=[],power_species=[]):
        self.initialize_class(file,L,As,ns,linear_species,power_species)
