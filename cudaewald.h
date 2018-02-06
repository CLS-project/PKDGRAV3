#ifndef CUDAEWALD_H
#define CUDAEWALD_H
#ifdef USE_CUDA
#include "cudautil.h"
#ifdef __cplusplus
extern "C" {
#endif
    int CUDA_queueEwald(void *cudaCtx,workParticle *work);
#ifdef __cplusplus
    }
#endif
#endif
#endif
