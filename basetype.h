#ifndef BASETYPE_H
#define BASETYPE_H
#include <stdint.h>
#include "ilp.h"
#include "ilc.h"

#ifdef HAVE_ALLOCA
#define stack_alloc(s) alloca(s)
#define stack_free(p) do {} while(0)
#else
#define stack_alloc(s) malloc(s)
#define stack_free(p) free(p)
#endif

#define IORDERBITS 41
#define IORDERMAX ((((uint64_t) 1)<<IORDERBITS)-1)

typedef struct particle {
    /*-----Base-Particle-Data----*/
    uint64_t iOrder     :  IORDERBITS;
    uint8_t  bMarked    :  1;
    uint8_t  uNewRung   :  6;
    uint8_t  uRung      :  6;
    uint8_t  bSrcActive :  1;
    uint8_t  bDstActive :  1;
    uint8_t  iClass     :  8;
    double r[3];
#ifdef SMOOTH_VARS
    /*-----Used-for-Smooth-------*/
    float fBallPRIVATE;
    float fDensityPRIVATE;
#endif
    } PARTICLE;

typedef float vel_t;

#define PP_CUDA_MEMORY_LIMIT (1024*1024)

typedef struct {
    float r[3];
    float a[3];
    float fSmooth2;
    float fDensity;
/*    float v[3];*/
/*    float fMass;*/
/*    float fSoft;*/
    } PINFOIN;

typedef struct {
    float a[3];
    float fPot;
    float dirsum, normsum;
    float rhopmax;
    } PINFOOUT;

typedef union {
#if defined(USE_SIMD)
    float *f;
#ifndef __CUDACC__
    v_sf *p;
#endif
#else
    double *f;
#endif
    } ewaldFloat;

typedef struct {
    ewaldFloat hx,hy,hz;
    ewaldFloat hCfac,hSfac;
    } EwaldTable;
struct EwaldVariables {
    double r[3]; /* Center of mass of the box */
    MOMC mom; /* moment of the box */
    double fEwCut2,fInner2,alpha,ialpha,alpha2,k1,ka,Lbox;
    double Q4xx,Q4xy,Q4xz,Q4yy,Q4yz,Q4zz,Q4,Q3x,Q3y,Q3z,Q2;
    int nMaxEwhLoop;
    int nEwLoopInner, nEwhLoop;
    int nReps,nEwReps;
    };
/*
** Accumulates the work for a set of particles
*/
typedef struct {
    PARTICLE **pPart;
    PINFOIN *pInfoIn;
    PINFOOUT *pInfoOut;
    float dRhoFac;
    int nP;
    int nRefs;
    void *ctx;
    int bGravStep;
#ifdef USE_CUDA
    void *cudaCtx;
#endif
    } workParticle;

/*
** One tile of PP interactions
*/
typedef struct {
    PINFOOUT *pInfoOut;
    ILP ilp;
    ILPTILE tile;
    workParticle *work;
    int i;
    } workPP;

typedef struct {
    PINFOOUT *pInfoOut;
    ILC ilc;
    ILCTILE tile;
    workParticle *work;
    int i;
    } workPC;

/* Careful! For compute <3.0, 65535 is the limit */
#define MAX_EWALD_PARTICLES 16384
typedef struct {
    workParticle **ppWorkPart;
    int *piWorkPart;
    void * pkd;
    int nP;
    } workEwald;


#endif
