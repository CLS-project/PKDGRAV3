#ifndef CLUTIL_H
#define CLUTIL_H
#ifdef __cplusplus
#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif
#include "basetype.h"

#define CL_PP_MAX_BUFFERED 128


typedef struct cl_wq_node {
    /* We can put this on different types of queues */
    struct cl_wq_node *next;
    void *ctx;
    int (*initFcn)(void *,void *,void *);
    int (*checkFcn)(void *,void *);
    cl_mem memOutCPU, memInGPU, memOutGPU, memInCPU;
    void *pHostBufToGPU, *pHostBufFromGPU;
    void *pClBufIn;
    void *pClBufOut;
    double startTime;
    cl_command_queue clQueue;
    cl_event clEvent;
    workParticle *ppWP[CL_PP_MAX_BUFFERED];
    int ppNI[CL_PP_MAX_BUFFERED];
    int ppSizeIn; // Number of bytes consumed in the buffer
    int ppSizeOut; // Number of bytes consumed in the buffer
    int ppnBlocks;
    int ppnBuffered;
    int bGravStep;
    union {
	struct {
	    size_t nBufferIn;
	    size_t nBufferOut;
	    int nGrid;
	    } pppc;
	};
    } CLwqNode;

typedef struct clContext {
    cl_context clContext;
    cl_device_id clDeviceId;
    cl_mem LxEwald,LyEwald,LzEwald,bhEwald, ewEwald, hxEwald, hyEwald, hzEwald, hCfac, hSfac;
    cl_program programEwald;
    } *CLCONTEXT;

typedef struct cl_ctx {
    int nCores, iCore;
    struct EwaldVariables *ewIn;
    EwaldTable *ewt;
    char hostname[256];
    int nWorkQueueSize, nWorkQueueBusy;
    int inClBufSize, outClBufSize;

    CLwqNode *wqCuda;
    CLwqNode *wqFree;
    CLCONTEXT context;

    cl_kernel kernelEwald;
    cl_command_queue queueEwald;
    cl_event eventEwald;

    } *CLCTX;

cl_program CL_compile(CLCTX cl, const char *src);
#else
void *CL_create_context();
void *CL_initialize(void *vctx,int nCores,int iCore);
void CL_SetQueueSize(void *vcl,int clSize, int inClBufSize, int outClBufSiz);
#endif
#endif
