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
//#include <gsl_odeiv.h>
//#include <gsl_sf_hyperg.h>

#ifdef CRAY_T3D
#include "hyperlib.h"
#endif
#include "cosmo.h"

#ifdef USE_GSL_COSMO
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
    //assert( csm->val.dOmegaDE == 0.0); // WARNING: no dark energy in this equation
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

double RK4_f(CSM csm, double lna, double G){
    double a = exp(lna);
    return G/csmExp2Hub(csm, a);
}

double RK4_g(CSM csm, double lna, double D, double G){
    double a = exp(lna);
    double inva = 1./a;
    return -2.0 * G + 1.5 * csm->val.dOmega0 * csm->val.dHubble0*csm->val.dHubble0 * inva*inva*inva * D/csmExp2Hub(csm, a);
}

/* THIS PART REQUIRES AN UPDATE OF THE GSL LIBRARY
 
int gsl_func(double a, const double y[], double f[], CSM csm){
    double inva = 1./a;
    double H = csmExp2Hub(csm,a);
    double invH = 1./H;
    f[0] = y[1]*invH;
    f[1] = -2 * y[1] + 1.5 * csm->val.dOmega0 * csm->val.dHubble0*csm->val.dHubble0 * inva*inva*inva * invH * y[0];
    return GSL_SUCCESS;
}

int gsl_jacobian(double a, const double y[], double *dfdy, double dfdt[], CSM csm){
    double H = csmExp2Hub(csm,a);
    double invH = 1./H;
    gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, 2, 2);
    gsl_matrix * m = &dfdy_mat.matrix;
    gsl_matrix_set (m, 0, 0, 0.0);
    gsl_matrix_set (m, 0, 1, invH);
    gsl_matrix_set (m, 1, 0, 1.5 * csm->val.dOmega0 * csm->val.dHubble0*csm->val.dHubble0 * inva*inva*inva * invH);
    gsl_matrix_set (m, 1, 1, 2.0);

    dfdt[0] = 0.0;
    dfdt[1] = 0.0;

   return GSL_SUCCESS;
}
*/

void MyRK4(CSM csm, double a, float *D1, float *f1){
   
    printf("Cosmology: Omega0 = %.3f, OmegaRad = %.5g, OmegaDE = %.3f, w0 = %.3f, wa = %.3f\n", csm->val.dOmega0, csm->val.dOmegaRad,csm->val.dOmegaDE,csm->val.w0,csm->val.wa);

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
	printf("\nSolving ODE using RK4...\n");

	double k0, k1, k2, k3;
        double l0, l1, l2, l3;

        FILE *f = fopen("GrowthFactors_RK4.txt", "w");

        if (f==NULL) {
	    printf("Error opening growth factor data file!\n");
            exit(1);
        }

        fprintf(f, "#a, z,D(z)\n");
         
        int i; // running loop variable
        for(i=0;i<=nSteps;i++){
	    ln_timesteps[i] = lna_init + i*stepwidth;
        fprintf(f, "%.15f, %.5f,%.20f\n", exp(ln_timesteps[i]),1.0/exp(ln_timesteps[i])-1.0, D[i]);
 
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
        
        fclose(f);

	*D1 = D[nSteps];
        *f1 = G[nSteps]/(csmExp2Hub(csm,a) * *D1); 

    return;

}

void csmComoveGrowth(CSM csm, double a, float *GrowthFactor, float *GrowthRate) {
    double dOmegaRad = csm->val.dOmegaRad;
    float GrowthFactorLCDM;
    double RadiativeCorrection_D;
    double RadiativeCorrection_f;

    /* 
    ** START: In this section we compute D_LCDM and hence intentionally set dOmegaRad = 0 because D_LCDM by definition has no radiation.
    */
    csm->val.dOmegaRad = 0;
    
    double eta = csmExp2Hub(csm, a);
    double result = ComoveGrowthFactorIntegral(csm,a);
 
    MyRK4(csm, a, &GrowthFactorLCDM, GrowthRate);
    
    printf("\nFrom new RK4 solver:\n===================\nGrowthFactor_LCDM = %.15f\n", GrowthFactorLCDM);

    /*
    ** END
    */

    /* 
    ** For the computation of radiative correction term we set dOmegaRad back to its original value.
    */
    csm->val.dOmegaRad = dOmegaRad;
    RadiativeCorrection_D = 2.0/3.0*csmRadMatEquivalence(csm);
    printf("Radiative Correction = %.10f\n", RadiativeCorrection_D);
 
    *GrowthFactor = GrowthFactorLCDM + RadiativeCorrection_D;

    RadiativeCorrection_f = 2.0/3.0*csmRadMatEquivalence(csm)/ *GrowthFactor; 

    *GrowthRate = *GrowthRate + RadiativeCorrection_f;   

    return;
    }

double csmComoveGrowthFactor(CSM csm,double a) {
    double dOmegaRad = csm->val.dOmegaRad;
    double GrowthFactorLCDM;
    double RadiativeCorrection;
    
    /*
    ** START: In this section we compute D_LCDM and hence intentionally set dOmegaRad = 0 because D_LCDM by definition has no radiation.
    */
    csm->val.dOmegaRad = 0;
    
    double eta = csmExp2Hub(csm, a);
    double result = ComoveGrowthFactorIntegral(csm,a);
    
    GrowthFactorLCDM = csm->val.dHubble0*csm->val.dHubble0*2.5*csm->val.dOmega0*eta*result;
    
    printf("\nFrom old PKDRGAV:\n=================\nGrowthFactor_LCDM = %.15f\n", GrowthFactorLCDM);

    /*
    ** END
    */
    
    /*
    ** For the computation of radiative correction term we set dOmegaRad back to its original value.
    */
    csm->val.dOmegaRad = dOmegaRad;
    RadiativeCorrection = 2.0/3.0*csmRadMatEquivalence(csm);
    
    printf("Radiative Correction = %.10f\n", RadiativeCorrection);

    return GrowthFactorLCDM + RadiativeCorrection;
}

double csmComoveGrowthRate(CSM csm,double a) {
    double dOmegaRad = csm->val.dOmegaRad;
    double D1;
    double GrowthRateLCDM;
    double RadiativeCorrection;

    /*
    ** START: In this section we compute D_LCDM and hence intentionally set dOmegaRad = 0 because D_LCDM by definition has no radiation.
    */
    csm->val.dOmegaRad = 0;
    GrowthRateLCDM = csmExp2HubRate(csm,a) + a*csmComoveGrowthInt(csm,a)/ComoveGrowthFactorIntegral(csm,a);
    /*
    ** END
    */
    
    /*
    ** For the computation of radiative correction term we set dOmegaRad back to its original value.
    */
    csm->val.dOmegaRad = dOmegaRad;
    RadiativeCorrection = 2.0/3.0*csmRadMatEquivalence(csm)/csmComoveGrowthFactor(csm,a);
    return (1 - RadiativeCorrection) * GrowthRateLCDM;
    }

void WriteGF(CSM csm, double a) {

    double a_init = 1e-12;
    int nSteps = 1000;
    double stepwidth = (a-a_init)/nSteps;
    
    FILE *f = fopen("GrowthFactors_PKDGRAV.txt", "w");

    if (f==NULL) {
       printf("Error opening growth factor data file!\n");
       exit(1);
    }

    fprintf(f, "#a,z,D(z)\n");

    int counter;
    double a_step;
    for (counter = 0; counter <= nSteps; counter++){
        a_step = a_init + counter*stepwidth;
	fprintf(f, "%.15f, %.5f, %.20f\n", a_step, 1.0/a_step -1.0, csmComoveGrowthFactor(csm, a_step));	
    }

    fclose(f);

    return; 
}
