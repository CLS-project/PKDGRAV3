#ifndef HYDRO_H
#define HYDRO_H
#include "pkd.h"
#include "smooth/smoothfcn.h"
#include <vector>

#define XX 0
#define YY 3
#define ZZ 5
#define XY 1
#define XZ 2
#define YZ 4

template<typename T>
int sign(T v) {return (v>0) - (v<0);}

typedef double my_real;

/* -----------------
 * MAIN FUNCTIONS AND CLASSES
 * -----------------
 */

/* Density loop */
struct hydroDensityPack {
    blitz::TinyVector<double,3> position;
    float fBall;
    uint8_t iClass;
    bool bMarked;
};

void hydroDensity(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf);
void hydroDensity_node(PKD pkd, SMF *smf, Bound bnd_node, const std::vector<PARTICLE *> &sinks,
                       NN *nnList, int nCnt);
void hydroDensityFinal(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf);
void packHydroDensity(void *vpkd,void *dst,const void *src);
void unpackHydroDensity(void *vpkd,void *dst,const void *src);


/* Gradient loop */
struct hydroGradientsPack {
    blitz::TinyVector<double,3> position;
    blitz::TinyVector<double,3> velocity;
    double P;
    float fBall;
    float fDensity;
    uint8_t iClass;
    bool bMarked;
};

void hydroGradients(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf);
void packHydroGradients(void *vpkd,void *dst,const void *src);
void unpackHydroGradients(void *vpkd,void *dst,const void *src);


/* Flux loop */
struct hydroFluxesPack {
    blitz::TinyVector<double,3> position;
    blitz::TinyVector<double,3> velocity;
    blitz::TinyVector<double,6> B;
    blitz::TinyVector<myreal,3> gradRho, gradVx, gradVy, gradVz, gradP;
    myreal lastUpdateTime;
    blitz::TinyVector<myreal,3> lastAcc;
    double omega;
    double P;
    float fBall;
    float fDensity;
    uint8_t uRung;
    uint8_t iClass;
    bool bMarked;
};

struct hydroFluxesFlush {
    myreal Frho;
    blitz::TinyVector<myreal,3> Fmom;
    myreal Fene;
#ifndef USE_MFM
    blitz::TinyVector<double,3> drDotFrho;
#endif
    blitz::TinyVector<double,3> mom;
    double E;
    double Uint;
    float fMass;
};

void hydroRiemann(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf);
void hydroRiemann_vec(PARTICLE *p,float fBall,int nSmooth, int nBuff,
                      my_real *restrict input_buffer,
                      my_real *restrict output_buffer,
                      SMF *smf);
void packHydroFluxes(void *vpkd,void *dst,const void *src);
void unpackHydroFluxes(void *vpkd,void *dst,const void *src);
void initHydroFluxes(void *vpkd,void *vp);
void initHydroFluxesCached(void *vpkd,void *vp);
void flushHydroFluxes(void *vpkd,void *dst,const void *src);
void combHydroFluxes(void *vpkd,void *p1,const void *p2);
void hydroFluxFillBuffer(my_real *input_buffer, PARTICLE *q, int i, int nBuff,
                         double dr2, blitz::TinyVector<double,3> dr, SMF *);
void hydroFluxUpdateFromBuffer(my_real *output_buffer, my_real *input_buffer,
                               PARTICLE *p, PARTICLE *q, int i, int nBuff, SMF *);
void hydroFluxGetNvars(int *in, int *out);

void pkdResetFluxes(PKD pkd,double dTime,double dDelta,double,double);


/* Time step loop */
struct hydroStepPack {
    blitz::TinyVector<double,3> position;
    blitz::TinyVector<double,3> velocity;
    float c;
    float fBall;
    uint8_t uRung;
    uint8_t uNewRung;
    uint8_t iClass;
    bool bMarked;
};

struct hydroStepFlush {
    uint8_t uNewRung;
    uint8_t uWake;
};

void hydroStep(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf);
void packHydroStep(void *vpkd,void *dst,const void *src);
void unpackHydroStep(void *vpkd,void *dst,const void *src);
void initHydroStep(void *vpkd,void *dst);
void flushHydroStep(void *vpkd,void *dst,const void *src);
void combHydroStep(void *vpkd,void *dst,const void *src);
void pkdWakeParticles(PKD pkd,int iRoot,double dTime,double dDelta);


/* Source terms */
void hydroSourceGravity(PKD pkd, particleStore::ParticleReference &p, SPHFIELDS *psph,
                        double pDelta, blitz::TinyVector<double,3> &pa, double dScaleFactor,
                        int bComove);
void hydroSourceExpansion(PKD pkd, particleStore::ParticleReference &p, SPHFIELDS *psph,
                          double pDelta, double dScaleFactor, double dHubble,
                          int bComove, double dConstGamma);
void hydroSyncEnergies(PKD pkd, particleStore::ParticleReference &p, SPHFIELDS *psph,
                       const blitz::TinyVector<double,3> &pa, double dConstGamma);
void hydroSetPrimitives(PKD pkd, particleStore::ParticleReference &p, SPHFIELDS *psph,
                        double dTuFac, double dConstGamma);
void hydroSetLastVars(PKD pkd, particleStore::ParticleReference &p, SPHFIELDS *psph,
                      const blitz::TinyVector<double,3> &pa, double dScaleFactor,
                      double dTime, double dDelta, double dConstGamma);

/* -----------------
 * HELPERS
 * -----------------
 */
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
void BarthJespersenLimiter(double *limVar, const blitz::TinyVector<double,3> &gradVar,
                           double var_max, double var_min,
                           const blitz::TinyVector<double,3> &dr);
void ConditionedBarthJespersenLimiter(double *limVar, const blitz::TinyVector<myreal,3> &gradVar,
                                      double var_max, double var_min,
                                      const blitz::TinyVector<double,3> &dr,
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

        phi_min = std::min(Lstate, Rstate);
        phi_max = std::max(Lstate, Rstate);

        if (sign(phi_min - d1) == sign(phi_min) ) {
            phi_m = phi_min - d1;
        }
        else {
            phi_m = phi_min/(1. + d1/fabs(phi_min));
        }

        if (sign(phi_max + d1) == sign(phi_max) ) {
            phi_p = phi_max + d1;
        }
        else {
            phi_p = phi_max/(1. + d1/fabs(phi_max));
        }

        if (Lstate < Rstate) {
            *Lstate_face = std::max(phi_m, std::min(phi_mean+d2, *Lstate_face));
            *Rstate_face = std::min(phi_p, std::max(phi_mean-d2, *Rstate_face));
        }
        else {
            *Rstate_face = std::max(phi_m, std::min(phi_mean+d2, *Rstate_face));
            *Lstate_face = std::min(phi_p, std::max(phi_mean-d2, *Lstate_face));
        }

    }


}
void compute_Ustar(double rho_K, double S_K, double v_K,
                   double p_K, double h_K, double S_s,
                   double *rho_sK, double *rhov_sK, double *e_sK);

#endif
