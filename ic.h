#ifndef IC_H
#define IC_H
#ifdef MDL_FFTW
typedef union {
    double *r;
    fftw_complex *k;
    } gridptr;

void pkdGenerateIC(PKD pkd,MDLFFT fft,int iBegYr,int iEndYr,int iBegZk,int iEndZk,gridptr dic[]);
#endif
#endif
