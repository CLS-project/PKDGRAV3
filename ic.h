#ifndef IC_H
#define IC_H
#ifdef MDL_FFTW
typedef union {
    FFTW3(real) *r;
    FFTW3(complex) *k;
    } gridptr;

typedef struct {
    FFTW3(real) x,y,z;
    } gridpos;

typedef struct {
    FFTW3(real) x,y,z;
    FFTW3(real) vx,vy,vz;
    } gridpsc;

void pkdGenerateIC(PKD pkd,int iSeed,double dBoxSize,double dOmegaMatter,double dOmegaLambda,double a,
    MDLFFT fft,int iBegYr,int iEndYr,int iBegZk,int iEndZk,gridptr dic[],gridpos *pos);
#endif
#endif
