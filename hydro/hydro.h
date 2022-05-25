#ifndef HYDRO_H
#define HYDRO_H
#include "pkd.h"
#include "smooth/smoothfcn.h"


#define XX 0
#define YY 3
#define ZZ 5
#define XY 1
#define XZ 2
#define YZ 4

typedef double my_real;

/* -----------------
 * MAIN FUNCTIONS
 * -----------------
 */
#ifdef __cplusplus
extern "C" {
#endif

/* Density loop */
void hydroDensity(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf);
void hydroDensity_node(PKD pkd, SMF *smf, BND bnd_node, PARTICLE **sinks, NN *nnList,
                       int nCnt_own, int nCnt);
void hydroDensityFinal(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf);

/* Gradient loop */
void hydroGradients(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf);

/* Flux loop */
void initHydroFluxes(void *vpkd, void *vp);
void initHydroFluxesCached(void *vpkd, void *vp);
void hydroRiemann(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf);
void hydroRiemann_vec(PARTICLE *p,float fBall,int nSmooth,
                      my_real **restrict input_buffer,
                      my_real **restrict output_buffer,
                      SMF *smf);
void pkdResetFluxes(PKD pkd,double dTime,double dDelta,double,double);

void combThirdHydroLoop(void *vpkd, void *p1,const void *p2);
void hydroFluxFillBuffer(my_real **buffer, PARTICLE *q, int i,
                         double dr2, double dx, double dy, double dz, SMF *);
void hydroFluxUpdateFromBuffer(my_real **out_buffer, my_real **in_buffer,
                               PARTICLE *p, PARTICLE *q, int i, SMF *);
void hydroFluxGetNvars(int *in, int *out);

/* Time step loop */
void initHydroStep(void *vpkd, void *vp);
void hydroStep(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf);
void combHydroStep(void *vpkd, void *p1,const void *p2);
void pkdWakeParticles(PKD pkd,int iRoot, double dTime, double dDelta);


/* Source terms */
void hydroSourceGravity(PKD pkd, PARTICLE *p, SPHFIELDS *psph,
                        double pDelta, double *pa, double dScaleFactor,
                        int bComove);
void hydroSourceExpansion(PKD pkd, PARTICLE *p, SPHFIELDS *psph,
                          double pDelta, double dScaleFactor, double dHubble,
                          int bComove, double dConstGamma);
void hydroSyncEnergies(PKD pkd, PARTICLE *p, SPHFIELDS *psph, double pa[3],
                       double dConstGamma);

void hydroSetPrimitives(PKD pkd, PARTICLE *p, SPHFIELDS *psph, double dTuFac, double dConstGamma);

void hydroSetLastVars(PKD pkd, PARTICLE *p, SPHFIELDS *psph, double *pa,
                      double dScaleFactor, double dTime, double dDelta,
                      double dConstGamma);

/* -----------------
 * HELPERS
 * -----------------
 */
#define SIGN(x) (((x) > 0) ? 1 : (((x) < 0) ? -1 : 0) )
#define MIN(X, Y)  ((X) < (Y) ? (X) : (Y))
#define MAX(X, Y)  ((X) > (Y) ? (X) : (Y))
void inverseMatrix(double *E, double *B);
double conditionNumber(double *E, double *B);
inline double cubicSplineKernel(double r, double h) {
    double q;
    q = r/h;
    if (q<1.0) {
        return M_1_PI/(h*h*h)*( 1. - 1.5*q*q*(1.-0.5*q) );
    }
    else if (q<2.0) {
        return 0.25*M_1_PI/(h*h*h)*(2.-q)*(2.-q)*(2.-q);
    }
    else {
        return 0.0;
    }
}
void BarthJespersenLimiter(double *limVar, double *gradVar,
                           double var_max, double var_min,
                           double dx, double dy, double dz);
void ConditionedBarthJespersenLimiter(double *limVar, myreal *gradVar,
                                      double var_max, double var_min,
                                      double dx, double dy, double dz,
                                      double Ncrit, double Ncond);
#define psi1 0.5
#define psi2 0.25
#pragma omp declare simd
inline void genericPairwiseLimiter(double Lstate, double Rstate,
                                   double *Lstate_face, double *Rstate_face) {
#ifdef DEBUG_FLUX_NOLIMITER
    return;
#endif
    double phi_max, phi_min, d1, d2, phi_mean, phi_p, phi_m;

    if (Lstate == Rstate) {
        *Lstate_face = Lstate;
        *Rstate_face = Rstate;
    }
    else {

        d1 = psi1*fabs(Lstate - Rstate);
        d2 = psi2*fabs(Lstate - Rstate);

        phi_mean = 0.5*(Lstate+Rstate);

        phi_min = MIN(Lstate, Rstate);
        phi_max = MAX(Lstate, Rstate);

        if (SIGN(phi_min - d1) == SIGN(phi_min) ) {
            phi_m = phi_min - d1;
        }
        else {
            phi_m = phi_min/(1. + d1/fabs(phi_min));
        }

        if (SIGN(phi_max + d1) == SIGN(phi_max) ) {
            phi_p = phi_max + d1;
        }
        else {
            phi_p = phi_max/(1. + d1/fabs(phi_max));
        }

        if (Lstate < Rstate) {
            *Lstate_face = MAX(phi_m, MIN(phi_mean+d2, *Lstate_face));
            *Rstate_face = MIN(phi_p, MAX(phi_mean-d2, *Rstate_face));
        }
        else {
            *Rstate_face = MAX(phi_m, MIN(phi_mean+d2, *Rstate_face));
            *Lstate_face = MIN(phi_p, MAX(phi_mean-d2, *Lstate_face));
        }

    }


}
void compute_Ustar(double rho_K, double S_K, double v_K,
                   double p_K, double h_K, double S_s,
                   double *rho_sK, double *rhov_sK, double *e_sK);

#ifdef __cplusplus
}
#endif
#endif
