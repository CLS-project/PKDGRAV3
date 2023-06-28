# distutils: language = c++
# cython: always_allow_keywords=True

import cython
from math import sqrt,pi
import numpy as np

class Cosmology:
    """
    This is the base class for cosmology. It is not meant to be constructed directly,
    but rather through one of the subclasses :class:`SimpleCosmology` or :class:`ClassCosmology`.
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

    @property
    def running(self):
        """
        :return: running of the spectral index :math:`dn_s/d\\ln k`
        :rtype: number
        """
        return self._csm.val.dRunning
    
    @property
    def pivot(self):
        """
        :return: pivot scale for the power spectrum :math:`k_p`
        :rtype: number
        """
        return self._csm.val.dPivot
    
    @property
    def comoving(self):
        """
        :return: whether to use comoving or physical units
        :rtype: bool
        """
        return self._csm.val.bComove
    
    @property
    def sigma8(self):
        """
        :return: sigma8 :math:`\\sigma_8`
        :rtype: number
        """
        return self._csm.val.dSigma8
    
    @property
    def As(self):
        """
        :return: amplitude of the power spectrum :math:`A_s`
        :rtype: number
        """
        return self._csm.val.dNormalization
    
    @property
    def ns(self):
        """
        :return: spectral index :math:`n_s`
        :rtype: number
        """
        return self._csm.val.dSpectral
    
    @property
    def h(self):
        """
        :return: hubble parameter :math:`h`
        :rtype: number
        """
        return self._csm.val.h
    
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
    :param number k: array of wavenumbers :math:`k`
    :param number Tk: array of transfer functions :math:`T(k)`
    :param str transfer_file: file name of the transfer function (input to numpy.loadtxt)
    :param tuple usecols: columns to use in the transfer function file (passed to numpy.loadtxt)
    """
    def __init__(self,*,omega_matter,omega_lambda=None,omega_baryon=0.0,omega_radiation=0.0,omega_dark_energy=0.0,
                        w0=-1.0,wa=0.0,running=0.0,pivot=0.05,comoving=True,sigma8=0.0,As=0.0,ns=0.0,h=1.0,
                        k=None,Tk=None,transfer_file=None,usecols=(0,1)):
        if transfer_file is not None:
            if k is not None or Tk is not None:
                raise ValueError("Cannot specify both transfer_file and k,Tk")
            if len(usecols) != 2:
                raise ValueError("usecols must be a tuple of length 2")
            k,Tk = np.loadtxt(transfer_file,usecols=usecols,unpack=True)
        self._k = k
        self._Tk = Tk
        self.initialize_simple(omega_matter,omega_lambda if omega_lambda is not None else 1.0 - omega_matter - omega_radiation - omega_dark_energy,
                               omega_baryon,omega_radiation,omega_dark_energy,w0,wa,running,pivot,comoving,sigma8,As,ns,h,k,Tk)
    def __repr__(self):
        return (f"SimpleCosmology("
                f"omega_matter={self.omega_matter},omega_lambda={self.omega_lambda},omega_baryon={self.omega_baryon},"
                f"omega_radiation={self.omega_radiation},omega_dark_energy={self.omega_dark_energy},w0={self.w0},wa={self.wa},"
                f"running={self.running},pivot={self.pivot},comoving={self.comoving},"
                f"sigma8={self.sigma8},As={self.As},ns={self.ns},h={self.h},"
                f"k={repr(self._k)},Tk={repr(self._Tk)})")

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
        self._file = file
        self._L = L
        self._linear_species = linear_species
        self._power_species = power_species
        self.initialize_class(file,L,As,ns,linear_species,power_species)

    def __repr__(self):
        return (f"ClassCosmology("
                f"file={self._file},L={self._L},As={self.As},ns={self.ns},"
                f"linear_species={repr(self._linear_species)},power_species={repr(self._power_species)})")

        # return cosmo.csmComoveKickFac(self._csm, time, delta);
    def rho_bar_matter(self,a):
        """
        Mean density :math:`\\bar{\\rho}` of matter.

        :param number a: expansion factor :math:`a`
        :return: mean matter density :math:`\\bar{\\rho}(a)`
        :rtype: number
        """
        return cosmo.csmRhoBar_m(self._csm,a)

    def rho_bar_linear(self,a):
        """
        Mean density :math:`\\bar{\\rho}` of linear species.

        :param number a: expansion factor :math:`a`
        :return: mean density :math:`\\bar{\\rho}(a)` of linear species
        :rtype: number
        """
        return cosmo.csmRhoBar_lin(self._csm,a)

    def rho_bar_power(self,a):
        """
        Mean density :math:`\\bar{\\rho}` of linear species to include in the power spectrum.

        :param number a: expansion factor :math:`a`
        :return: mean density :math:`\\bar{\\rho}(a)` of linear species to include in the power spectrum
        :rtype: number
        """
        return cosmo.csmRhoBar_pk(self._csm,a)


    def delta_matter(self,a,k):
        """
        The :math:`\\delta` perturbation of matter species.

        :param number a: expansion factor :math:`a`
        :param number k: wavenumber :math:`k`
        :return: delta perturbations :math:`\\delta(a,k)` of matter species
        :rtype: number
        """
        return cosmo.csmDelta_m(self._csm,a,k)

    def theta_matter(self,a,k):
        """
        The :math:`\\theta` perturbation of matter species.

        :param number a: expansion factor :math:`a`
        :param number k: wavenumber :math:`k`
        :return: delta perturbations :math:`\\theta(a,k)` of matter species
        :rtype: number
        """
        return cosmo.csmTheta_m(self._csm,a,k)

    def delta_linear(self,a,k):
        """
        The :math:`\\delta` perturbation of linear species.

        :param number a: expansion factor :math:`a`
        :param number k: wavenumber :math:`k`
        :return: delta perturbations :math:`\\delta(a,k)` of linear species
        :rtype: number
        """
        return cosmo.csmDelta_lin(self._csm,a,k)

    def delta_power(self,a,k):
        """
        The :math:`\\delta` perturbation of linear species to include in the power spectrum.

        :param number a: expansion factor :math:`a`
        :param number k: wavenumber :math:`k`
        :return: delta perturbations :math:`\\delta(a,k)` of linear species to include in the power spectrum
        :rtype: number
        """
        return cosmo.csmDelta_pk(self._csm,a,k)

    def delta_rho_linear(self,a,a_next,k):
        """
        Calculates :math:`\\delta\\bar{\\rho}` of linear species.

        :param number a: expansion factor :math:`a`
        :param number a_next: next expansion factor :math:`a'`
        :param number k: wavenumber :math:`k`
        :return: :math:`\\delta(a,k)\\bar{\\rho}(a)` of linear species
        :rtype: number
        """
        return cosmo.csmDeltaRho_lin(self._csm,a,a_next,k)

    def delta_rho_power(self,a,k):
        """
        Calculates :math:`\\delta\\bar{\\rho}` for the species included in the power spectrum.

        :param number a: expansion factor :math:`a`
        :param number k: wavenumber :math:`k`
        :return: :math:`\\delta(a,k)\\bar{\\rho}(a)`
        :rtype: number
        """
        return cosmo.csmDeltaRho_pk(self._csm,a,k)

    def zeta(self,k):
        """
        :param number k: wavenumber :math:`k`
        :return: zeta :math:`\\zeta(k)`
        :rtype: number
        """
        return cosmo.csmZeta(self._csm,k)
