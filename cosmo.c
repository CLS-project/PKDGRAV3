#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <ctype.h>
#include <getopt.h>
#include <assert.h>

#ifdef CRAY_T3D
#include "hyperlib.h"
#endif
#include "cosmo.h"

//This version of the code requires GSL
#ifndef USE_GSL_COSMO
#error USE_GSL_COSMO must be defined!
#else
#define LIMIT 1000
#endif

#define OPT_z   'z'
#define OPT_rad 'r'
#define OPT_Om0 'm'
#define OPT_w0 'w'
//#define OPT_wa 'wa'
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
    csm->W = gsl_integration_workspace_alloc(LIMIT);
    *pcsm = csm;
    }

void csmFinish(CSM csm) {
    gsl_integration_workspace_free(csm->W);
    free(csm);
    }

#define EPSCOSMO 1e-7

/*
 ** by MK: Computes the scale factor a at radiation-matter equivalence.
 */
double csmRadMatEquivalence(CSM csm){
    return csm->val.dOmegaRad/csm->val.dOmega0;
}

double csmTime2Hub(CSM csm,double dTime) {
    double a = csmTime2Exp(csm,dTime);

    assert(a > 0.0);
    return csmExp2Hub(csm, a);
    }

static double Exp2Time_integrand(double ak, void *params) {
    CSM csm = (CSM)params;

    double dExp = pow(ak,2.0/3.0);
    assert (dExp > 0.0);

    return 2.0/(3.0*ak*csmExp2Hub(csm,dExp));
    }

static double Exp2TimeIntegrate(CSM csm,double dExp) {
    gsl_function F;
    F.function = &Exp2Time_integrand;
    F.params = csm;
    double result,error;
    gsl_integration_qag(&F, 0.0, pow(dExp, 1.5),
	0.0, EPSCOSMO, LIMIT, GSL_INTEG_GAUSS61, csm->W, &result, &error);
    printf("a=%g,\t result of Exp2TimeIntegrate = %g\n", dExp, result);
    return result;
    }

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
    if (csm->val.dOmegaDE == 0.0 && csm->val.dLambda == 0.0 && csm->val.dOmegaRad == 0.0) {
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
	printf("Enter GSL integration...\n");
	double r = Exp2TimeIntegrate(csm,dExp);
	printf("The result of csmExp2Time is: %g\n",r);
	return r;
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

	f = dTime - Exp2TimeIntegrate(csm,a);
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

	    f = dTime - Exp2TimeIntegrate(csm,a);
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
    else if (csm->val.dOmegaDE == 0.0 && csm->val.dLambda == 0.0 && csm->val.dOmegaRad == 0.0) {
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
	gsl_function F;
	F.function = &ComoveDrift_integrand;
	F.params = csm;
	double result,error;
	gsl_integration_qag(&F, 
            1.0/csmTime2Exp(csm, dTime), 1.0/csmTime2Exp(csm, dTime + dDelta),
            0.0, EPSCOSMO, LIMIT, GSL_INTEG_GAUSS61, csm->W, &result, &error);
        return result;
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
    else if (csm->val.dOmegaDE == 0.0 && csm->val.dLambda == 0.0 && csm->val.dOmegaRad == 0.0) {
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
	gsl_function F;
	F.function = &ComoveKick_integrand;
	F.params = csm;
	double result,error;
	gsl_integration_qag(&F, 
	    1.0/csmTime2Exp(csm, dTime), 1.0/csmTime2Exp(csm, dTime + dDelta),
	    0.0, EPSCOSMO, LIMIT, GSL_INTEG_GAUSS61, csm->W, &result, &error);
	return result;
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

double RK4_f(CSM csm, double lna, double G){
    double a = exp(lna);
    return G/csmExp2Hub(csm, a);
}


double RK4_g(CSM csm, double lna, double D, double G){
    double a = exp(lna);
    double inva = 1./a;
    return -2.0 * G + 1.5 * csm->val.dOmega0 * csm->val.dHubble0*csm->val.dHubble0 * inva*inva*inva * D/csmExp2Hub(csm, a);
}

void MyRK4(CSM csm, double a, float *D1, float *f1){
    /*
    ** Variable declarations & initializations
    */
    double a_init, lna_init = log(1e-12); // ln(a)=-12 ==> a = e^(-12) ~ 0 

    int nSteps = 1000;
    double stepwidth = (log(a)- lna_init)/nSteps;

    // NOTICE: Storing the following quantities into data structures is by far not optimal (we actually never need the old values after the update).
    double ln_timesteps[nSteps+1];
    float D[nSteps+1]; // Growth factor D(a)
    float G[nSteps+1]; // G(a) = dD(a)/dln(a) *H  ==>  Growth rate: f(a) = G/(H*D) 

    /* 
    ** Set boundary conditions
    */
    a_init = exp(lna_init);
    D[0] = a_init;
    G[0] = csmExp2Hub(csm, a_init)*a_init;

    //Classical RK4 Solver - This is the best of all solvers available in this code!
    double k0, k1, k2, k3;
    double l0, l1, l2, l3;

    int i; // running loop variable
    for(i=0;i<=nSteps;i++){
	ln_timesteps[i] = lna_init + i*stepwidth;
   
    	//RK4 step 1
       	k0 = stepwidth * RK4_f(csm, ln_timesteps[i], G[i]);
 	l0 = stepwidth * RK4_g(csm, ln_timesteps[i], D[i], G[i]);

	//RK4 step 2
	k1 = stepwidth * RK4_f(csm, ln_timesteps[i] + stepwidth/2.0, G[i] + l0/2.0); 
	l1 = stepwidth * RK4_g(csm, ln_timesteps[i] + stepwidth/2.0, D[i] + k0/2.0, G[i] + l0/2.0);
    
       	//RK4 step 3
	k2 = stepwidth * RK4_f(csm, ln_timesteps[i] + stepwidth/2.0, G[i] + l1/2.0);
	l2 = stepwidth * RK4_g(csm, ln_timesteps[i] + stepwidth/2.0, D[i] + k1/2.0, G[i] + l1/2.0);

	//RK4 step 4
	k3 = stepwidth * RK4_f(csm, ln_timesteps[i] + stepwidth, G[i] + l2);
	l3 = stepwidth * RK4_g(csm, ln_timesteps[i] + stepwidth, D[i] + k2, G[i] + l2);

	//Update
	D[i+1] = D[i] + (k0 + 2*k1 + 2*k2 + k3)/6.0;
	G[i+1] = G[i] + (l0 + 2*l1 + 2*l2 + l3)/6.0; 
    }
        
    *D1 = D[nSteps];
    *f1 = G[nSteps]/(csmExp2Hub(csm,a) * *D1); 

    return;
}
