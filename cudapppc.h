#ifndef CUDAPPPC_H
#define CUDAPPPC_H
#ifdef USE_CUDA
#include "cudautil.h"
#include "ilp.h"
#include "ilc.h"
#ifdef __cplusplus
extern "C" {
#endif
    int CUDA_queuePP(void *cudaCtx,workParticle *wp, ILPTILE tile, int bGravStep);
    int CUDA_queuePC(void *cudaCtx,workParticle *wp, ILCTILE tile, int bGravStep);
#ifdef __cplusplus
    }
#endif
#endif
#endif
