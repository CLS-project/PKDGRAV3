#ifndef IC_H
#define IC_H
#include <complex.h>

typedef struct {
    float r[3];
    float v[3];
    } basicParticle;

#ifdef MDL_FFTW

typedef union {
    FFTW3(real) *r;
//    FFTW3(complex) *k;
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
    double dOmega0,double dLambda0,double dSigma8,double dSpectral,double h,
    double a,int nTf, double *tk, double *tf,
    double *noiseMean, double *noiseCSQ);
#endif
#endif
