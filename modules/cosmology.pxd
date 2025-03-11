from libc.stdlib cimport malloc, free
cimport cosmo
cimport numpy as np
from math import sqrt,pi

cdef extern from "Python.h":
    const char* PyUnicode_AsUTF8(object unicode)

cdef class Cosmology:
    cdef cosmo.csmContext* _csm

    cdef inline initialize(self):
        cosmo.csmInitialize(&self._csm)

    cdef inline finish(self):
        if self._csm is not NULL:
            cosmo.csmFinish(self._csm)

    cdef inline cosmo.csmContext * csm(self):
        return self._csm

    cdef inline _comoving_growth(self,a):
        cdef double D1LPT, D2LPT, D3aLPT, D3bLPT, D3cLPT, f1LPT, f2LPT, f3aLPT, f3bLPT, f3cLPT
        cosmo.csmComoveGrowth(self._csm, a,
                              &D1LPT, &D2LPT, &D3aLPT, &D3bLPT, &D3cLPT,
                              &f1LPT, &f2LPT, &f3aLPT, &f3bLPT, &f3cLPT)
        return (D1LPT, D2LPT, D3aLPT, D3bLPT, D3cLPT, f1LPT, f2LPT, f3aLPT, f3bLPT, f3cLPT)




cdef class SimpleCosmology(Cosmology):
    cdef np.ndarray _k
    cdef np.ndarray _Tk

    cdef inline initialize_simple(self,omega_matter,omega_lambda,omega_baryon,omega_radiation,omega_dark_energy,
                        w0,wa,running,pivot,comoving,sigma8,As,ns,h,k,Tk):
        self._csm.val.classData.bClass = 0
        self._csm.val.bComove = comoving
        self._csm.val.dHubble0 = sqrt(8.0/3.0 * pi)
        self._csm.val.dOmega0 = omega_matter
        self._csm.val.dLambda = omega_lambda
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
    cdef object _file
    cdef object _linear_species
    cdef object _power_species
    cdef double _L

    cdef inline const char** parse_species(self, list species, int n):
        cdef:
            int i
            const char **c_species = <const char**> malloc(sizeof(char*) * n)
    
        if c_species is NULL:
            raise MemoryError()
            
        for i in range(n):
            c_species[i] = PyUnicode_AsUTF8(species[i])
        
        return c_species

    cdef inline initialize_class(self,file,L,As,ns,linear_species,power_species):
        
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
