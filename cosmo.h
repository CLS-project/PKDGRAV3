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

#ifndef COSMO_HINCLUDED
#define COSMO_HINCLUDED

#define USE_GSL_COSMO
#ifdef USE_GSL_COSMO
    #include <gsl/gsl_integration.h>
    #include <gsl/gsl_spline.h>
    #include <gsl/gsl_interp2d.h>
    #include <gsl/gsl_spline2d.h>
#endif

/*
** Nested strucst storing the CLASS data.
** They are nested in the following way:
** - classDataStruct (classData)
**   - classDataBackgroundStruct (background)
**   - classDataPerturbationsStruct (perturbations)
*/
#define CLASS_BACKGROUND_SIZE 10000
#define CLASS_PERTURBATIONS_A_SIZE 1024
#define CLASS_PERTURBATIONS_K_SIZE 256
struct classDataBackgroundStruct {
    size_t size;
    double a      [CLASS_BACKGROUND_SIZE];
    double t      [CLASS_BACKGROUND_SIZE];
    double H      [CLASS_BACKGROUND_SIZE];
    double rho_m  [CLASS_BACKGROUND_SIZE];
    double rho_lin[CLASS_BACKGROUND_SIZE];
    double rho_pk [CLASS_BACKGROUND_SIZE];
};
struct classDataPerturbationsStruct {
    size_t size_a;
    size_t size_k;
    double a[CLASS_PERTURBATIONS_A_SIZE];
    double k[CLASS_PERTURBATIONS_K_SIZE];
    double delta_m  [CLASS_PERTURBATIONS_A_SIZE*CLASS_PERTURBATIONS_K_SIZE];
    double theta_m  [CLASS_PERTURBATIONS_A_SIZE*CLASS_PERTURBATIONS_K_SIZE];
    double delta_lin[CLASS_PERTURBATIONS_A_SIZE*CLASS_PERTURBATIONS_K_SIZE];
    double delta_pk [CLASS_PERTURBATIONS_A_SIZE*CLASS_PERTURBATIONS_K_SIZE];
};
struct classDataStruct {
    int bClass;
    int nLinear; /* Number of linear species */
    int nPower; /* Number of linear species for measuring P(k) */
    struct classDataBackgroundStruct background;
    struct classDataPerturbationsStruct perturbations;
};

/* Nested strucst storing the CLASS GSL objects.
** They are nested in the following way:
** - classGslStruct (classGsl)
**   - classGslBackgroundStruct (background)
**   - classGslPerturbationsStruct (perturbations)
*/
struct classGslBackgroundStruct {
    gsl_interp_accel *logExp2logHub_acc;
    gsl_spline       *logExp2logHub_spline;
    gsl_interp_accel *logTime2logHub_acc;
    gsl_spline       *logTime2logHub_spline;
    gsl_interp_accel *logExp2logTime_acc;
    gsl_spline       *logExp2logTime_spline;
    gsl_interp_accel *logTime2logExp_acc;
    gsl_spline       *logTime2logExp_spline;
    gsl_interp_accel *logExp2logRho_m_acc;
    gsl_spline       *logExp2logRho_m_spline;
    gsl_interp_accel *logExp2logRho_lin_acc;
    gsl_spline       *logExp2logRho_lin_spline;
    gsl_interp_accel *logExp2logRho_pk_acc;
    gsl_spline       *logExp2logRho_pk_spline;
};
struct classGslPerturbationsStruct {
    gsl_interp_accel *logk2delta_m_acc;
    gsl_interp_accel *loga2delta_m_acc;
    gsl_spline2d     *logkloga2delta_m_spline;
    gsl_interp_accel *logk2theta_m_acc;
    gsl_interp_accel *loga2theta_m_acc;
    gsl_spline2d     *logkloga2theta_m_spline;
    gsl_interp_accel *logk2delta_lin_acc;
    gsl_interp_accel *loga2delta_lin_acc;
    gsl_spline2d     *logkloga2delta_lin_spline;
    gsl_interp_accel *logk2delta_pk_acc;
    gsl_interp_accel *loga2delta_pk_acc;
    gsl_spline2d     *logkloga2delta_pk_spline;
};
struct classGslStruct {
    int initialized;
    struct classGslBackgroundStruct background;
    struct classGslPerturbationsStruct perturbations;
};

struct csmVariables {
    int bComove;
    double dHubble0;
    double dOmega0;
    double dLambda;
    double dOmegaRad;
    double dOmegab;
    double dOmegaDE;
    double w0;
    double wa;
    double dSigma8;
    double dNormalization;  /* either sigma8 or normalization must be non-zero */
    double dSpectral;
    double dRunning;
    double dPivot;
    double h;
    struct classDataStruct classData;
};


typedef struct csmContext {
    struct csmVariables val;

#ifdef USE_GSL_COSMO
    gsl_integration_workspace *W;
#endif
    struct classGslStruct classGsl;
} *CSM;
#ifdef __cplusplus
extern "C" {
#endif
void csmClassRead(CSM csm, const char *achFilename, double dBoxSize, double h,
                  int nLinear, const char **aLinear, int nPower, const char **aPower);
void csmClassGslInitialize(CSM csm);
double csmRhoBar_m    (CSM csm, double a);
double csmRhoBar_lin  (CSM csm, double a);
double csmRhoBar_pk   (CSM csm, double a);
double csmDelta_m     (CSM csm, double a,                double k);
double csmTheta_m     (CSM csm, double a,                double k);
double csmDelta_lin   (CSM csm, double a,                double k);
double csmDelta_pk    (CSM csm, double a,                double k);
double csmDeltaRho_lin(CSM csm, double a, double a_next, double k);
double csmDeltaRho_pk (CSM csm, double a,                double k);
double csmZeta        (CSM csm,                          double k);

void csmInitialize(CSM *pcsm);
void csmFinish(CSM csm);
double csmRadMatEquivalence(CSM csm);

static inline double csmExp2Hub(CSM csm, double dExp) {
    if (csm->val.classData.bClass) {
        if (dExp > csm->val.classData.background.a[csm->val.classData.background.size - 1]) {
            /* dExp is in the future; do linear extrapolation */
            return csm->val.classData.background.H[csm->val.classData.background.size - 1]
                   + (
                       csm->val.classData.background.H[csm->val.classData.background.size - 1]
                       - csm->val.classData.background.H[csm->val.classData.background.size - 2]
                   )/(
                       csm->val.classData.background.a[csm->val.classData.background.size - 1]
                       - csm->val.classData.background.a[csm->val.classData.background.size - 2]
                   )*(dExp - csm->val.classData.background.a[csm->val.classData.background.size - 1]);
        }
        return exp(gsl_spline_eval(
                       csm->classGsl.background.logExp2logHub_spline,
                       log(dExp),
                       csm->classGsl.background.logExp2logHub_acc));
    }

    double dOmegaCurve = 1.0 - csm->val.dOmega0 - csm->val.dLambda - csm->val.dOmegaDE - csm->val.dOmegaRad;

    assert(dExp > 0.0);
    return csm->val.dHubble0
           *sqrt(csm->val.dOmega0*dExp
                 + dOmegaCurve*dExp*dExp
                 + csm->val.dOmegaRad
                 + csm->val.dOmegaDE*pow(dExp,1.0 - 3.0*(csm->val.w0 + csm->val.wa))*exp(-3.0*csm->val.wa*(1.0 - dExp))
                 + csm->val.dLambda*dExp*dExp*dExp*dExp)/(dExp*dExp);
}

double csmTime2Hub(CSM csm, double dTime);
double csmExp2Time(CSM csm, double dExp);
double csmTime2Exp(CSM csm, double dTime);
double csmComoveDriftInt(CSM csm, double dIExp);
double csmComoveKickInt(CSM csm, double dIExp);
double csmComoveDriftFac(CSM csm, double dTime, double dDelta);
double csmComoveKickFac(CSM csm, double dTime, double dDelta);
double csmComoveLookbackTime2Exp(CSM csm, double dComoveTime);
void csmComoveGrowth(CSM csm, double a, double *D1LPT, double *D2LPT, double *f1LPT, double *f2LPT);
#ifdef __cplusplus
}
#endif

/*
 ** returns the speed of light in simulation units, given
 ** the simulation length unit in h^-1 Mpc.
 */
static inline double dLightSpeedSim(double dMpcUnit) {
    /*
    ** Find the speed of light in simulation units.
    **
    ** c[Mpc/Gyr] = c[cm/s] * Julian Year[s] / pc[cm] * 1000
    ** c_sim = c[Mpc/Gyr] * (x Gyrs/ 1 sim time) * ( 1 sim length/Boxsize (Mpc))
    ** x = 1/sqrt(4.498*h*h*2.776e-4)
    */
    /*return(8676.85/dMpcUnit);*/

    /*
    ** Doug's version:
    **
    ** Cosmological coordinates
    ** G     = 4.30172e-9 Mpc/M. (km/s)^2
    ** rho_c = 3 H^2 / (8 pi G)
    ** c     = 299792.458 km/s
    **
    ** c_sim = c[km/s] * sqrt(Lbox / (G * rho_c * Lbox^3))
    **       = c[km/s] * sqrt(8 pi / (3 H^2 Lbox^2) )
    **       = c[km/s] * sqrt(8 pi / 3) / Lbox / H
    **       = c[km/s] * sqrt(8 pi / 3) / Lbox / h / 100
    ** dMpcUnit given in Mpc/h gives:
    **       = 299792.458 * sqrt(8 pi / 3) / 100 / dMpcUnit
    **       = 8677.2079486362706 / dMpcUnit
    */
    return 8677.2079486362706 / dMpcUnit;
}

#endif
