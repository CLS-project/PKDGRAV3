# distutils: language = c++
# cython: always_allow_keywords=True

from libc.stdlib cimport malloc, free
cimport cosmo
from math import sqrt,pi

cdef extern from "Python.h":
    const char* PyUnicode_AsUTF8(object unicode)

cdef class Cosmology:
    cdef cosmo.csmContext* _csm

    def __cinit__(self):
        cosmo.csmInitialize(&self._csm)
        pass
    def __dealloc__(self):
        if self._csm is not NULL:
            cosmo.csmFinish(self._csm)
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
        cdef double D1LPT, D2LPT, f1LPT, f2LPT
        cosmo.csmComoveGrowth(self._csm, a, &D1LPT, &D2LPT, &f1LPT, &f2LPT)
        return (D1LPT, D2LPT, f1LPT, f2LPT)


cdef class SimpleCosmology(Cosmology):
    def __init__(self,*,omega_matter,omega_lambda=None,omega_baryon=0.0,omega_radiation=0.0,omega_dark_energy=0.0,
                        w0=-1.0,wa=0.0,running=0.0,pivot=0.05,comoving=True,sigma8=0.0,As=0.0,ns=0.0,h=1.0):
        self._csm.val.classData.bClass = 0
        self._csm.val.bComove = comoving
        self._csm.val.dHubble0 = sqrt(8.0/3.0 * pi)
        self._csm.val.dOmega0 = omega_matter
        self._csm.val.dLambda = omega_lambda if omega_lambda is not None else 1.0 - omega_matter - omega_radiation - omega_dark_energy
        self._csm.val.dOmegaRad = omega_radiation
        self._csm.val.dOmegab = omega_baryon
        self._csm.val.dOmegaDE = omega_dark_energy
        self._csm.val.w0 = w0
        self._csm.val.wa = wa
        self._csm.val.dSigma8 = sigma8
        self._csm.val.dNormalization = As
        self._csm.val.dSpectral = ns
        self._csm.val.dRunning = running
        self._csm.val.dPivot = pivot
        self._csm.val.h = h

cdef class ClassCosmology(Cosmology):
    cdef const char** parse_species(self, list species, int n):
        cdef:
            int i
            const char **c_species = <const char**> malloc(sizeof(char*) * n)
    
        if c_species is NULL:
            raise MemoryError()
            
        for i in range(n):
            c_species[i] = PyUnicode_AsUTF8(species[i])
        
        return c_species

    def __init__(self,file,L=1.0,As=0.0,ns=0.0,linear_species=[],power_species=[]):
        cdef:
            int n_linear = len(linear_species)
            int n_power = len(power_species)
            const char *c_file
            const char **c_linear_species
            const char **c_power_species

        c_file = PyUnicode_AsUTF8(file)
        
        try:
            c_linear_species = self.parse_species(linear_species, n_linear)
            c_power_species = self.parse_species(power_species, n_power)

            self._csm.val.classData.bClass = 1
            self._csm.val.bComove = 1
            self._csm.val.dHubble0 = sqrt(8.0/3.0 * pi)
            self._csm.val.dNormalization = As
            self._csm.val.dSpectral = ns
            cosmo.csmClassRead(self._csm, c_file, L, 0.0, n_linear,c_linear_species,n_power,c_power_species)
            cosmo.csmClassGslInitialize(self._csm);

        finally:
            free(c_linear_species)
            free(c_power_species)