# distutils: language = c++
# cython: always_allow_keywords=True

import cython
from math import sqrt,pi

class Cosmology:
    """
    This is the base class for cosmology. It is not meant to be used directly.
    """
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
        """
        :return: matter density :math:`\\Omega_m`
        :rtype: number
        """
        return self._csm.val.dOmega0

    @property
    def omega_baryon(self):
        """
        :return: baryon density :math:`\\Omega_b`
        :rtype: number
        """
        return self._csm.val.dOmegab

    @property
    def omega_lambda(self):
        """
        :return: dark energy density :math:`\\Omega_\\Lambda`
        :rtype: number
        """
        return self._csm.val.dLambda

    @property
    def omega_radiation(self):
        """
        :return: radiation density :math:`\\Omega_\\gamma`
        :rtype: number
        """
        return self._csm.val.dOmegaRad

    @property
    def omega_dark_energy(self):
        """
        :return: dark energy density :math:`\\Omega_{DE}`
        :rtype: number
        """
        return self._csm.val.dOmegaDE

    @property
    def w0(self):
        """
        :return: dark energy equation of state parameter :math:`w_0`
        :rtype: number
        """
        return self._csm.val.w0

    @property
    def wa(self):
        """
        :return: dark energy equation of state parameter :math:`w_a`
        :rtype: number
        """
        return self._csm.val.wa

    def radiation_matter_equivalence(self):
        """
        :return: radiation-matter equivalence :math:`a_{eq}`
        :rtype: number
        """
        return cosmo.csmRadMatEquivalence(self._csm)
        
    def expansion_to_hubble(self,a):
        """
        :param number a: expansion factor :math:`a`
        :return: hubble parameter :math:`H(a)`
        :rtype: number
        """
        return cosmo.csmExp2Hub(self._csm,a)

    def time_to_hubble(self,time):
        """
        :param number time: time :math:`t`
        :return: hubble parameter :math:`H(t)`
        :rtype: number
        """
        return cosmo.csmTime2Hub(self._csm,time)

    def expansion_to_time(self,a):
        """
        :param number a: expansion factor :math:`a`
        :return: time :math:`t(a)`
        :rtype: number
        """
        return cosmo.csmExp2Time(self._csm,a)

    def time_to_expansion(self,time):
        """
        :param number time: time :math:`t`
        :return: expansion factor :math:`a(t)`
        :rtype: number
        """
        return cosmo.csmTime2Exp(self._csm,time)

    def comoving_drift_factor(self,time,delta):
        """
        :param number time: time :math:`t`
        :param number delta: delta time :math:`\\Delta t`
        :return: comoving drift factor
        :rtype: number
        """
        return cosmo.csmComoveDriftFac(self._csm, time, delta);

    def comoving_kick_factor(self,time,delta):
        """
        :param number time: time :math:`t`
        :param number delta: delta time :math:`\\Delta t`
        :return: comoving kick factor
        :rtype: number
        """
        return cosmo.csmComoveKickFac(self._csm, time, delta);

    def comoving_growth(self,a):
        """
        :param number a: expansion factor :math:`a`
        :return: comoving growth
        :rtype: number
        """
        return self._comoving_growth(a)

class SimpleCosmology(Cosmology):
    """
    This is a simple cosmology class. Use this by setting
    
    :param number omega_matter: matter density :math:`\\Omega_m`
    :param number omega_lambda: dark energy density :math:`\\Omega_\\Lambda`
    :param number omega_baryon: baryon density :math:`\\Omega_b`
    :param number omega_radiation: radiation density :math:`\\Omega_\\gamma`
    :param number omega_dark_energy: dark energy density :math:`\\Omega_{DE}`
    :param number w0: dark energy equation of state parameter :math:`w_0`
    :param number wa: dark energy equation of state parameter :math:`w_a`
    :param number running: running of the spectral index :math:`dn_s/d\\ln k`
    :param number pivot: pivot scale for the power spectrum :math:`k_p`
    :param bool comoving: whether to use comoving or physical units
    :param number sigma8: sigma8 :math:`\\sigma_8`
    :param number As: amplitude of the power spectrum :math:`A_s`
    :param number ns: spectral index :math:`n_s`
    :param number h: hubble parameter :math:`h`
    """
    def __init__(self,*,omega_matter,omega_lambda=None,omega_baryon=0.0,omega_radiation=0.0,omega_dark_energy=0.0,
                        w0=-1.0,wa=0.0,running=0.0,pivot=0.05,comoving=True,sigma8=0.0,As=0.0,ns=0.0,h=1.0):
        self.initialize_simple(omega_matter,omega_lambda if omega_lambda is not None else 1.0 - omega_matter - omega_radiation - omega_dark_energy,
                               omega_baryon,omega_radiation,omega_dark_energy,w0,wa,running,pivot,comoving,sigma8,As,ns,h)

class ClassCosmology(Cosmology):
    """
    This is a cosmology class that uses CLASS. Use this by setting
    
    :param str file: File name of the CLASS/Concept HDF5 file
    :param number L: box size in Mpc/h
    :param number As: amplitude of the power spectrum :math:`A_s`
    :param number ns: spectral index :math:`n_s`
    :param list linear_species: list of linear species
    :param list power_species: list of power species
    """
    def __init__(self,file,L=1.0,As=0.0,ns=0.0,linear_species=[],power_species=[]):
        self.initialize_class(file,L,As,ns,linear_species,power_species)
