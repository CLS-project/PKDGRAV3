#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#ifdef CRAY_T3D
#include "hyperlib.h"
#endif

#include "cosmo.h"

#ifdef USE_GSL_COSMO
#define LIMIT 1000
#endif

/*
 * Cosmological module for PKDGRAV.
 * N.B.  This code is being shared with skid and the I.C. generator.
 */

void csmInitialize(CSM *pcsm) {
    CSM csm;

    csm = (CSM) malloc(sizeof(struct csmContext));
    assert(csm != NULL);

    csm->val.dHubble0 = 0.0;
    csm->val.dOmega0 = 0.0;
    csm->val.dLambda = 0.0;
    csm->val.dOmegaDE = 0.0;
    csm->val.w0 = 0.0;
    csm->val.wa = 0.0;
    csm->val.dOmegaRad = 0.0;
    csm->val.dOmegab = 0.0;
    csm->val.bComove = 0;
#ifdef USE_GSL_COSMO
    csm->W = gsl_integration_workspace_alloc(LIMIT);
#endif
    *pcsm = csm;
    }

void csmFinish(CSM csm) {
#ifdef USE_GSL_COSMO
    gsl_integration_workspace_free(csm->W);
#endif
    free(csm);
    }

#define EPSCOSMO 1e-7

#ifndef USE_GSL_COSMO
double dRombergO(void *CTX, double (*func)(void *, double), double a,
		 double b, double eps);
#endif

/*
 * ** by MK: Computes the scale factor a at radiation-matter equivalence.
 * */
double csmRadMatEquivalence(CSM csm){
    return csm->val.dOmegaRad/csm->val.dOmega0;
}

/*
 ** The cosmological equation of state is entirely determined here.  We
 ** will derive all other quantities from this function.
 */

double csmExp2Hub(CSM csm, double dExp) {
    double dOmegaCurve = 1.0 - csm->val.dOmega0 -
			 csm->val.dLambda - csm->val.dOmegaDE - csm->val.dOmegaRad;

    assert(dExp > 0.0);
    return csm->val.dHubble0
	   *sqrt(csm->val.dOmega0*dExp
		 + dOmegaCurve*dExp*dExp
		 + csm->val.dOmegaRad
		 + csm->val.dOmegaDE*pow(dExp,1.0 - 3.0*(csm->val.w0 + csm->val.wa))*exp(-3.0*csm->val.wa*(1.0 - dExp))
		 + csm->val.dLambda*dExp*dExp*dExp*dExp)/(dExp*dExp);
    }

/*
 * ** This is a/H(a)*dH(a)/da
 * */
double csmExp2HubRate(CSM csm, double dExp) {
    double dOmegaCurve = 1.0 - csm->val.dOmega0 -
			 csm->val.dLambda - csm->val.dOmegaDE - csm->val.dOmegaRad;
    assert(dExp > 0.0);
    double ia = 1.0 / dExp;
    assert( csm->val.dOmegaDE == 0.0); // WARNING: no dark energy in this equation
    return ( -1.5 * csm->val.dOmega0 * ia*ia*ia - dOmegaCurve * ia*ia - 2.0*csm->val.dOmegaRad * ia*ia*ia*ia )
	 / (        csm->val.dOmega0 * ia*ia*ia + dOmegaCurve * ia*ia +     csm->val.dOmegaRad * ia*ia*ia*ia + csm->val.dLambda );
    }


double csmTime2Hub(CSM csm,double dTime) {
    double a = csmTime2Exp(csm,dTime);

    assert(a > 0.0);
    return csmExp2Hub(csm, a);
    }

double csmCosmoTint(CSM csm, double dY) {
    double dExp = pow(dY, 2.0/3.0);

    assert(dExp > 0.0);
    return 2.0/(3.0*dY*csmExp2Hub(csm, dExp));
    }

#ifdef USE_GSL_COSMO
static double Exp2Time_integrand(double ak, void * params) {
    return csmCosmoTint(params,ak);
    }
static double Exp2TimeIntegrate(CSM csm,double dExp) {
    gsl_function F;
    F.function = &Exp2Time_integrand;
    F.params = csm;
    double result,error;
    gsl_integration_qag(&F, 0.0, pow(dExp, 1.5),
	0.0, EPSCOSMO, LIMIT, GSL_INTEG_GAUSS61, csm->W, &result, &error);
    return result;
    }
#endif

double csmExp2Time(CSM csm,double dExp) {
    double dOmega0 = csm->val.dOmega0;
    double dHubble0 = csm->val.dHubble0;
    double a0,A,B,eta;

    if (!csm->val.bComove) {
	/*
	 ** Invalid call!
	 */
	assert(0);
	}
    if (csm->val.dLambda == 0.0 && csm->val.dOmegaRad == 0.0) {
	if (dOmega0 == 1.0) {
	    assert(dHubble0 > 0.0);
	    if (dExp == 0.0) return(0.0);
	    return(2.0/(3.0*dHubble0)*pow(dExp,1.5));
	    }
	else if (dOmega0 > 1.0) {
	    assert(dHubble0 >= 0.0);
	    if (dHubble0 == 0.0) {
		B = 1.0/sqrt(dOmega0);
		eta = acos(1.0-dExp);
		return(B*(eta-sin(eta)));
		}
	    if (dExp == 0.0) return(0.0);
	    a0 = 1.0/dHubble0/sqrt(dOmega0-1.0);
	    A = 0.5*dOmega0/(dOmega0-1.0);
	    B = A*a0;
	    eta = acos(1.0-dExp/A);
	    return(B*(eta-sin(eta)));
	    }
	else if (dOmega0 > 0.0) {
	    assert(dHubble0 > 0.0);
	    if (dExp == 0.0) return(0.0);
	    a0 = 1.0/dHubble0/sqrt(1.0-dOmega0);
	    A = 0.5*dOmega0/(1.0-dOmega0);
	    B = A*a0;
	    eta = acosh(dExp/A+1.0);
	    return(B*(sinh(eta)-eta));
	    }
	else if (dOmega0 == 0.0) {
	    assert(dHubble0 > 0.0);
	    if (dExp == 0.0) return(0.0);
	    return(dExp/dHubble0);
	    }
	else {
	    /*
	     ** Bad value.
	     */
	    assert(0);
	    return(0.0);
	    }
	}
    else {
#ifdef USE_GSL_COSMO
	return Exp2TimeIntegrate(csm,dExp);
#else
	return dRombergO(csm, (double (*)(void *, double)) csmCosmoTint,
			 0.0, pow(dExp, 1.5), EPSCOSMO);
#endif
	}
    }

#define MAX_ITER 100

double csmTime2Exp(CSM csm,double dTime) {
    double al=0,ah=1,a0,a1=1,at,a;
    double th,f,f1,h,ho;
    int j;

    if (!csm->val.bComove) return(1.0);
    else {
	assert(dTime > 0);
	th = csmExp2Time(csm,ah);
	/*
	** Search for upper bracket if needed.
	*/
	while (dTime > th) {
	    a0 = a1;
	    a1 = ah;
	    ah = a1+a0;
	    th = csmExp2Time(csm,ah);
	    }
	a = 0.5*(al+ah);
	ho = ah-al;
	h = ho;
#ifdef USE_GSL_COSMO
	f = dTime - Exp2TimeIntegrate(csm,a);
#else
	f = dTime - dRombergO(csm, (double (*)(void *, double)) csmCosmoTint,0.0,pow(a,1.5),EPSCOSMO);
#endif
	f1 = 1/(a*csmExp2Hub(csm,a));
	for (j=0;j<MAX_ITER;++j) {
	    if (a+f/f1 < al || a+f/f1 > ah || fabs(2*f) > fabs(ho*f1)) {
		/*
		** Bisection Step.
		*/
		ho = h;
		h = 0.5*(ah-al);
		a = al+h;
		/*
				printf("bisect al:%.14g ah:%.14g a:%.14g\n",al,ah,a);
		*/
		if (a == al) return a;
		}
	    else {
		/*
		** Newton Step.
		*/
		ho = h;
		h = f/f1;
		at = a;
		a += h;
		/*
				printf("newton al:%.14g ah:%.14g a:%.14g\n",al,ah,a);
		*/
		if (a == at) return a;
		}
	    if (fabs(h) < EPSCOSMO) {
		/*
				printf("converged al:%.14g ah:%.14g a:%.14g t:%.14g == %.14g\n",
				       al,ah,a,dRombergO(csm, (double (*)(void *, double)) csmCosmoTint,0.0,pow(a,1.5),EPSCOSMO*1e-1),
				       dTime);
		*/
		return a;
		}
#ifdef USE_GSL_COSMO
	    f = dTime - Exp2TimeIntegrate(csm,a);
#else
	    f = dTime - dRombergO(csm, (double (*)(void *, double)) csmCosmoTint,0.0,pow(a,1.5),EPSCOSMO*1e-1);
#endif
	    f1 = 1/(a*csmExp2Hub(csm,a));
	    if (f < 0) ah = a;
	    else al = a;
	    }
	assert(0);
	}
    return 0.0; /* We never reach here, but this keeps the compiler happy */
    }


double csmComoveDriftInt(CSM csm, double dIExp) {
    return -dIExp/(csmExp2Hub(csm, 1.0/dIExp));
    }
static double ComoveDrift_integrand(double diExp, void * params) {
    return csmComoveDriftInt(params,diExp);
    }

/*
 ** Make the substitution y = 1/a to integrate da/(a^2*H(a))
 */
double csmComoveKickInt(CSM csm, double dIExp) {
    return -1.0/(csmExp2Hub(csm, 1.0/dIExp));
    }

static double ComoveKick_integrand(double diExp, void * params) {
    return csmComoveKickInt(params,diExp);
    }

/*
 ** This function integrates the time dependence of the "drift"-Hamiltonian.
 */
double csmComoveDriftFac(CSM csm,double dTime,double dDelta) {
    double dOmega0 = csm->val.dOmega0;
    double dHubble0 = csm->val.dHubble0;
    double a0,A,B,a1,a2,eta1,eta2;

    if (!csm->val.bComove) return(dDelta);
    else if (csm->val.dLambda == 0.0 && csm->val.dOmegaRad == 0.0) {
	a1 = csmTime2Exp(csm,dTime);
	a2 = csmTime2Exp(csm,dTime+dDelta);
	if (dOmega0 == 1.0) {
	    return((2.0/dHubble0)*(1.0/sqrt(a1) - 1.0/sqrt(a2)));
	    }
	else if (dOmega0 > 1.0) {
	    assert(dHubble0 >= 0.0);
	    if (dHubble0 == 0.0) {
		A = 1.0;
		B = 1.0/sqrt(dOmega0);
		}
	    else {
		a0 = 1.0/dHubble0/sqrt(dOmega0-1.0);
		A = 0.5*dOmega0/(dOmega0-1.0);
		B = A*a0;
		}
	    eta1 = acos(1.0-a1/A);
	    eta2 = acos(1.0-a2/A);
	    return(B/A/A*(1.0/tan(0.5*eta1) - 1.0/tan(0.5*eta2)));
	    }
	else if (dOmega0 > 0.0) {
	    assert(dHubble0 > 0.0);
	    a0 = 1.0/dHubble0/sqrt(1.0-dOmega0);
	    A = 0.5*dOmega0/(1.0-dOmega0);
	    B = A*a0;
	    eta1 = acosh(a1/A+1.0);
	    eta2 = acosh(a2/A+1.0);
	    return(B/A/A*(1.0/tanh(0.5*eta1) - 1.0/tanh(0.5*eta2)));
	    }
	else if (dOmega0 == 0.0) {
	    /*
	     ** YOU figure this one out!
	     */
	    assert(0);
	    return(0.0);
	    }
	else {
	    /*
	     ** Bad value?
	     */
	    assert(0);
	    return(0.0);
	    }
	}
    else {
#ifdef USE_GSL_COSMO
	gsl_function F;
	F.function = &ComoveDrift_integrand;
	F.params = csm;
	double result,error;
	gsl_integration_qag(&F, 
	    1.0/csmTime2Exp(csm, dTime), 1.0/csmTime2Exp(csm, dTime + dDelta),
	    0.0, EPSCOSMO, LIMIT, GSL_INTEG_GAUSS61, csm->W, &result, &error);
	return result;
#else
	return dRombergO(csm,
			 (double (*)(void *, double)) csmComoveDriftInt,
			 1.0/csmTime2Exp(csm, dTime),
			 1.0/csmTime2Exp(csm, dTime + dDelta), EPSCOSMO);
#endif
	}
    }


/*
 ** This function integrates the time dependence of the "kick"-Hamiltonian.
 */
double csmComoveKickFac(CSM csm,double dTime,double dDelta) {
    double dOmega0 = csm->val.dOmega0;
    double dHubble0 = csm->val.dHubble0;
    double a0,A,B,a1,a2,eta1,eta2;

    if (!csm->val.bComove) return(dDelta);
    else if (csm->val.dLambda == 0.0 && csm->val.dOmegaRad == 0.0) {
	a1 = csmTime2Exp(csm,dTime);
	a2 = csmTime2Exp(csm,dTime+dDelta);
	if (dOmega0 == 1.0) {
	    return((2.0/dHubble0)*(sqrt(a2) - sqrt(a1)));
	    }
	else if (dOmega0 > 1.0) {
	    assert(dHubble0 >= 0.0);
	    if (dHubble0 == 0.0) {
		A = 1.0;
		B = 1.0/sqrt(dOmega0);
		}
	    else {
		a0 = 1.0/dHubble0/sqrt(dOmega0-1.0);
		A = 0.5*dOmega0/(dOmega0-1.0);
		B = A*a0;
		}
	    eta1 = acos(1.0-a1/A);
	    eta2 = acos(1.0-a2/A);
	    return(B/A*(eta2 - eta1));
	    }
	else if (dOmega0 > 0.0) {
	    assert(dHubble0 > 0.0);
	    a0 = 1.0/dHubble0/sqrt(1.0-dOmega0);
	    A = 0.5*dOmega0/(1.0-dOmega0);
	    B = A*a0;
	    eta1 = acosh(a1/A+1.0);
	    eta2 = acosh(a2/A+1.0);
	    return(B/A*(eta2 - eta1));
	    }
	else if (dOmega0 == 0.0) {
	    /*
	     ** YOU figure this one out!
	     */
	    assert(0);
	    return(0.0);
	    }
	else {
	    /*
	     ** Bad value?
	     */
	    assert(0);
	    return(0.0);
	    }
	}
    else {
#ifdef USE_GSL_COSMO
	gsl_function F;
	F.function = &ComoveKick_integrand;
	F.params = csm;
	double result,error;
	gsl_integration_qag(&F, 
	    1.0/csmTime2Exp(csm, dTime), 1.0/csmTime2Exp(csm, dTime + dDelta),
	    0.0, EPSCOSMO, LIMIT, GSL_INTEG_GAUSS61, csm->W, &result, &error);
	return result;
#else
	return dRombergO(csm,
			 (double (*)(void *, double)) csmComoveKickInt,
			 1.0/csmTime2Exp(csm, dTime),
			 1.0/csmTime2Exp(csm, dTime + dDelta), EPSCOSMO);
#endif
	}
    }

double csmComoveLookbackTime2Exp(CSM csm,double dComoveTime) {
    if (!csm->val.bComove) return(1.0);
    else {
	double dExpOld = 0.0;
	double dT0 = csmExp2Time(csm, 1.0);
	double dTime = dT0 - dComoveTime;
	double dExpNew;
	int it = 0;

	if (dTime < EPSCOSMO) dTime = EPSCOSMO;
	dExpNew = csmTime2Exp(csm, dTime);
	/*
	 * Root find with Newton's method.
	 */
	do {
	    double dTimeNew = csmExp2Time(csm, dExpNew);
	    double f = dComoveTime
		       - csmComoveKickFac(csm, dTimeNew, dT0 - dTimeNew);
	    double fprime = -1.0/(dExpNew*dExpNew*csmExp2Hub(csm, dExpNew));
	    dExpOld = dExpNew;
	    dExpNew += f/fprime;
	    it++;
	    assert(it < 20);
	    }
	while (fabs(dExpNew - dExpOld)/dExpNew > EPSCOSMO);
	return dExpNew;
	}
    }

double csmComoveGrowthInt(CSM csm, double a) {
    if (a==0.0) return 0;
    else return pow(a*csmExp2Hub(csm,a),-3.0);
    }

static double ComoveGrowth_integrand(double a, void * params) {
    return csmComoveGrowthInt(params,a);
    }

static double ComoveGrowthFactorIntegral(CSM csm,double a) {
    double result;
#if defined(USE_GSL_COSMO)
    gsl_function F;
    F.function = &ComoveGrowth_integrand;
    F.params = csm;
    double error;
    gsl_integration_qag(&F, 
	0, a,0.0, 1e-12, LIMIT, GSL_INTEG_GAUSS61, csm->W, &result, &error);
#else
    result = dRombergO(csm,
	(double (*)(void *, double)) csmComoveGrowthInt,
	0, a, 1e-12);
#endif
    return result;
    }

double csmComoveGrowthFactor(CSM csm,double a) {
    // double a_equ = csm->val.dOmegaRad/csm->val.dOmega0; /* added by MK: computing the scale factor at radiation/matter equivalence */
    double dOmegaRad = csm->val.dOmegaRad;
    double GrowthFactorLCDM;
    double RadiativeCorrection;
    
    /* START */ /* In this section we compute D_LCDM and hence intentionally set dOmegaRad = 0 because D_LCDM by definition has no radiation. */
    csm->val.dOmegaRad = 0;

    double eta = csmExp2Hub(csm, a);
    double result = ComoveGrowthFactorIntegral(csm,a);
    
    GrowthFactorLCDM = csm->val.dHubble0*csm->val.dHubble0*2.5*csm->val.dOmega0*eta*result;
    /* END */
    /* For the radiative correction term we set dOmegaRad back to its original value */

    csm->val.dOmegaRad = dOmegaRad;
    RadiativeCorrection = 2.0/3.0*csmRadMatEquivalence(csm);    

    return GrowthFactorLCDM + RadiativeCorrection;
    }

double csmComoveGrowthRate(CSM csm,double a) {
    double dOmegaRad = csm->val.dOmegaRad;
    double D1;
    double GrowthRateLCDM;
    double RadiativeCorrection;

    /* START */ /* In this section we compute f1_LCDM and hence intentionally set dOmegaRad = 0 because f1_LCDM by definition has no radiation. */
    csm->val.dOmegaRad = 0;
    GrowthRateLCDM = csmExp2HubRate(csm,a) + a*csmComoveGrowthInt(csm,a)/ComoveGrowthFactorIntegral(csm,a);
    /* END */

    csm->val.dOmegaRad = dOmegaRad;
    RadiativeCorrection = 2.0/3.0*csmRadMatEquivalence(csm)/csmComoveGrowthFactor(csm,a);
    
    return (1 - RadiativeCorrection) * GrowthRateLCDM; 
    }


