#ifndef IC_H
#define IC_H

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
int pkdGenerateIC(PKD pkd,MDLFFT fft,int iSeed,int bFixed,float fPhase,int nGrid,int b2LPT,double dBoxSize,
    struct csmVariables *cosmo,double a,int nTf, double *tk, double *tf,
    double *noiseMean, double *noiseCSQ);
#endif
#endif
