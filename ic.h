#ifndef IC_H
#define IC_H
#include <complex.h>

typedef struct {
    float dr[3];
    float v[3];
    } basicParticle;

typedef struct {
    uint64_t ix : 21;
    uint64_t iy : 21;
    uint64_t iz : 21;
    float dr[3];
    float v[3];
    } expandParticle;

typedef union {
    basicParticle b;
    expandParticle e;
    } overlayedParticle;

#ifdef MDL_FFTW

typedef union {
    FFTW3(real) *r;
    float complex *k;
    } gridptr;

typedef struct {
    FFTW3(real) x,y,z;
    } gridpos;

typedef struct {
    FFTW3(real) x,y,z;
    FFTW3(real) vx,vy,vz;
    } gridpsc;

int pkdGenerateIC(PKD pkd,MDLFFT fft,int iSeed,int nGrid,int b2LPT,double dBoxSize,
    struct csmVariables *cosmo,double a,int nTf, double *tk, double *tf,
    double *noiseMean, double *noiseCSQ);
#endif
#endif
