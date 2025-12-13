#ifndef HYDRO_H
#define HYDRO_H
#include "core/simd.h"
#include "pkd.h"
#include "smooth/smoothfcn.h"
#include "riemann.h"
#include "limiters.h"
#include "solver.h"
#include "kernels.h"
#include <vector>



namespace meshless {
typedef double myreal;
}

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

void hydroDensity(particleStore::ParticleReference &p,float fBall,int nSmooth,NN *nnList,SMF *smf);
void hydroDensity_node(PKD pkd, SMF *smf, Bound bnd_node, const std::vector<PARTICLE *> &sinks,
                       NN *nnList, int nCnt);
void hydroDensityFinal(particleStore::ParticleReference &p,float fBall,int nSmooth,NN *nnList,SMF *smf);
class SmoothHelperHydroDensity : public SmoothHelperPacked<hydroDensityPack> {
    virtual void    pack(void *dst, const void *src) override final;
    virtual void  unpack(void *dst, const void *src, const void *key) override final;
    virtual void  init_shared(void *dst) override final {init_shared_gas_decorator(pkd,dst);}

public:
    SmoothHelperHydroDensity(pkdContext *pkd) : SmoothHelperPacked<hydroDensityPack>(pkd,false) {}
};

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

void hydroGradients(particleStore::ParticleReference &p,float fBall,int nSmooth,NN *nnList,SMF *smf);
class SmoothHelperHydroGradients : public SmoothHelperPacked<hydroGradientsPack> {
    virtual void    pack(void *dst, const void *src) override final;
    virtual void  unpack(void *dst, const void *src, const void *key) override final;
    virtual void  init_shared(void *dst) override final {init_shared_gas_decorator(pkd,dst);}
public:
    SmoothHelperHydroGradients(pkdContext *pkd) : SmoothHelperPacked<hydroGradientsPack>(pkd,false) {
        nCacheSize = pkd->decorators.DecoratorSize() + pkd->particles.ParticleSize();
    }
};

/* Flux loop */
struct hydroFluxesPack {
    blitz::TinyVector<double,3> position;
    blitz::TinyVector<double,3> velocity;
    blitz::TinyVector<double,6> B;
    blitz::TinyVector<meshless::myreal,3> gradRho, gradVx, gradVy, gradVz, gradP;
    meshless::myreal lastUpdateTime;
    blitz::TinyVector<meshless::myreal,3> lastAcc;
    double omega;
    double spheom;
    double P;
    float fBall;
    float fDensity;
    uint8_t uRung;
    uint8_t iClass;
    bool bMarked;
};

struct hydroFluxesFlush {
    meshless::myreal Frho;
    blitz::TinyVector<meshless::myreal,3> Fmom;
    meshless::myreal Fene;
    blitz::TinyVector<double,3> mom;
    double E;
    double Uint;
    // Used for the MFV method, ignored otherwise
    blitz::TinyVector<double,3> drDotFrho;
    float fMass;
};

void hydroRiemann(PARTICLE *p,float fBall,int nSmooth, int nBuff,
                  meshless::myreal *restrict input_buffer,
                  meshless::myreal *restrict output_buffer,
                  SMF *smf);
void hydroRiemann_wrapper(PARTICLE *p,float fBall,int nSmooth, int nBuff,
                          meshless::myreal *restrict input_buffer,
                          meshless::myreal *restrict output_buffer,
                          SMF *smf);
void hydroFluxFillBuffer(meshless::myreal *input_buffer, PARTICLE *q, int i, int nBuff,
                         double dr2, blitz::TinyVector<double,3> dr, SMF *);
void hydroFluxUpdateFromBuffer(meshless::myreal *output_buffer, meshless::myreal *input_buffer,
                               PARTICLE *p, PARTICLE *q, int i, int nBuff, SMF *);
void hydroFluxGetBufferInfo(int *in, int *out);

class SmoothHelperHydroFluxes : public SmoothHelperPacked<hydroFluxesPack,hydroFluxesFlush> {
    virtual void    pack(void *dst, const void *src) override final;
    virtual void  unpack(void *dst, const void *src, const void *key) override final;
    virtual void    init(void *dst) override final;
    virtual void   flush(void *dst, const void *src) override final;
    virtual void combine(void *dst, const void *src, const void *key) override final;
    virtual void  init_shared(void *dst) override final {init_shared_gas_decorator(pkd,dst);}
public:
    SmoothHelperHydroFluxes(pkdContext *pkd, bool bSymmetric=true) : SmoothHelperPacked<hydroFluxesPack,hydroFluxesFlush>(pkd,bSymmetric) {
        nCacheSize = pkd->decorators.DecoratorSize() + pkd->particles.ParticleSize();
    }
};

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

void hydroStep(particleStore::ParticleReference &p,float fBall,int nSmooth,NN *nnList,SMF *smf);
class SmoothHelperHydroStep : public SmoothHelperPacked<hydroStepPack,hydroStepFlush> {
    virtual void    pack(void *dst, const void *src) override final;
    virtual void  unpack(void *dst, const void *src, const void *key) override final;
    virtual void    init(void *dst) override final;
    virtual void   flush(void *dst, const void *src) override final;
    virtual void combine(void *dst, const void *src, const void *key) override final;
    virtual void  init_shared(void *dst) override final {init_shared_gas_decorator(pkd,dst);}
public:
    SmoothHelperHydroStep(pkdContext *pkd) : SmoothHelperPacked<hydroStepPack,hydroStepFlush>(pkd,true) {
        nCacheSize = pkd->decorators.DecoratorSize() + pkd->particles.ParticleSize();
    }
};


void pkdWakeParticles(PKD pkd,int iRoot,double dTime,double dDelta);

/* Source terms */
void hydroSourceGravity(PKD pkd, particleStore::ParticleReference &p, decorators::MFM *psph,
                        double pDelta, blitz::TinyVector<double,3> &pa, double dScaleFactor,
                        int bComove);
void hydroSourceExpansion(PKD pkd, particleStore::ParticleReference &p, decorators::MFM *psph,
                          double pDelta, double dScaleFactor, double dHubble,
                          int bComove, double dConstGamma);
void hydroSyncEnergies(PKD pkd, particleStore::ParticleReference &p, decorators::MFM *psph,
                       const blitz::TinyVector<double,3> &pa, double dConstGamma);
void hydroSetPrimitives(PKD pkd, particleStore::ParticleReference &p, decorators::MFM *psph,
                        double dTuFac, double dConstGamma);
void hydroSetLastVars(PKD pkd, particleStore::ParticleReference &p, decorators::MFM *psph,
                      const blitz::TinyVector<double,3> &pa, double dScaleFactor,
                      double dTime, double dDelta, double dConstGamma);
void hydroResetFluxes(decorators::MFM *psph);

/* -----------------
 * HELPERS
 * -----------------
 */
void inverseMatrix(double *E, double *B);
double conditionNumber(double *E, double *B);
void compute_Ustar(double rho_K, double S_K, double v_K,
                   double p_K, double h_K, double S_s,
                   double *rho_sK, double *rhov_sK, double *e_sK);


#endif
