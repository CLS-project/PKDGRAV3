/*  This file is part of PKDGRAV3 (http://www.pkdgrav.org/).
 *  Copyright (c) 2001-2018 Joachim Stadel & Douglas Potter
 *
 *  PKDGRAV3 is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  PKDGRAV3 is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with PKDGRAV3.  If not, see <http://www.gnu.org/licenses/>.
 */
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include "cosmo.h"

#define MAX_TF 500

typedef struct {
    gsl_interp_accel *acc;
    gsl_spline *spline;
    double *tk, *tf;
    double spectral;
    double normalization;
    int nTf;
    } powerParameters;

static double power(powerParameters *P,double k) {
    double T = gsl_spline_eval(P->spline,log(k),P->acc);
    return pow(k,P->spectral) * P->normalization * T * T;
    }

typedef struct {
    powerParameters *P;
    double r;
    } varianceParameters;

static double variance_integrand(double ak, void * params) {
    varianceParameters *vprm = params;
    double x, w;
    /* Window function for spherical tophat of given radius (e.g., 8 Mpc/h) */
    x = ak * vprm->r;
    w = 3.0*(sin(x)-x*cos(x))/(x*x*x);
    return power(vprm->P,ak)*ak*ak*w*w*4.0*M_PI;
    }

static double variance(powerParameters *P,double dRadius) {
    varianceParameters vprm;
    gsl_function F;
    double result, error;
    gsl_integration_workspace *W = gsl_integration_workspace_alloc (1000);
    vprm.P = P;
    vprm.r = dRadius; /* 8 Mpc/h for example */
    F.function = &variance_integrand;
    F.params = &vprm;
    gsl_integration_qag(&F, exp(P->tk[0]), exp(P->tk[P->nTf-1]),
	0.0, 1e-6, 1000, GSL_INTEG_GAUSS61, W, &result, &error);
    gsl_integration_workspace_free(W);
    return result;
    }

int main(int argc,char *argv[]) {
    CSM csm;
    int nTf;
    double k[MAX_TF];
    double tf[MAX_TF];
    char buffer[200];
    powerParameters P;
    double a;
    double D1_0, D2_0, f1_0, f2_0;
    double D1_a, D2_a, f1_a, f2_a;
    FILE *fp;

    if (argc<6) {
	fprintf(stderr, "Usage: %s tffile omega sigma8 spectral redshift\n",argv[0]);
	return 1;
	}
    csmInitialize(&csm);
    csm->val.dHubble0 = sqrt(8.0 * M_PI / 3.0);
    csm->val.dOmega0 = atof(argv[2]);
    csm->val.dLambda= 1.0 - csm->val.dOmega0;
    csm->val.dSigma8 = atof(argv[3]);
    csm->val.dSpectral = atof(argv[4]);
    a = 1.0 / (1.0 + atof(argv[5]));

    csmComoveGrowth(csm, 1.0, &D1_0, &D2_0, &f1_0, &f2_0);
    csmComoveGrowth(csm, a,   &D1_a, &D2_a, &f1_a, &f2_a);

    fp = fopen(argv[1],"r");
    if (fp == NULL) {
	perror(argv[1]);
	return 2;
	}

    nTf = 0;
    while(fgets(buffer,sizeof(buffer),fp)) {
	assert(nTf < MAX_TF);
	if (sscanf(buffer," %lg %lg\n",&k[nTf],&tf[nTf])==2) {
	    k[nTf] = log(k[nTf]);
	    ++nTf;
	    }
	}
    fclose(fp);


    P.normalization = 1.0;
    P.spectral = csm->val.dSpectral;
    P.nTf = nTf;
    P.tk = k;
    P.tf = tf;
    P.acc = gsl_interp_accel_alloc();
    P.spline = gsl_spline_alloc (gsl_interp_cspline, nTf);
    gsl_spline_init(P.spline, P.tk, P.tf, P.nTf);

    double dSigma8 = csm->val.dSigma8 * D1_a/D1_0;
    P.normalization *= dSigma8*dSigma8 / variance(&P,8.0);

    double twopi = 2.0 * 4.0 * atan(1.0);
    double twopi3 = pow(twopi,3.0);

    while(fgets(buffer,sizeof(buffer),stdin)) {
	double ak;
	if (sscanf(buffer,"%lg\n",&ak)!=1) {
	    fprintf(stderr,"ERROR\n");
	    return 3;
	    }
	double amp = power(&P,ak) * twopi3;
	printf("%g %g\n", ak, amp);
	}
    
    return 0;
    }
