#ifndef IC_H
#define IC_H
#ifdef MDL_FFTW
typedef union {
    double *r;
    fftw_complex *k;
    } gridptr;

typedef struct {
    double x,y,z;
    } gridpos;

typedef struct {
    double x,y,z;
    double vx,vy,vz;
    } gridpsc;

void pkdGenerateIC(PKD pkd,int iSeed,double dBoxSize,double dOmegaMatter,double dOmegaLambda,double a,
    MDLFFT fft,int iBegYr,int iEndYr,int iBegZk,int iEndZk,gridptr dic[],gridpos *pos);
#endif
#endif
