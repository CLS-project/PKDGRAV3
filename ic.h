#ifndef IC_H
#define IC_H

#ifdef __cplusplus
#include <complex>
typedef std::complex<float> COMPLEX;
#define REAL(x) std::real(x)
#define IMAG(x) std::imag(x)
static const std::complex<float> I = {0,1};
#else
#include <complex.h>
#define COMPLEX float complex
#define REAL(x) creal(x)
#define IMAG(x) cimag(x)
#endif

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

typedef struct {
    FFTW3(real) x,y,z;
    } gridpos;

typedef struct {
    FFTW3(real) x,y,z;
    FFTW3(real) vx,vy,vz;
    } gridpsc;

#ifdef __cplusplus
extern "C"
#endif
int pkdGenerateIC(PKD pkd,MDLFFT fft,int iSeed,int nGrid,int b2LPT,double dBoxSize,
    struct csmVariables *cosmo,double a,int nTf, double *tk, double *tf,
    double *noiseMean, double *noiseCSQ);
#endif
#endif
