#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <unistd.h>
#include <cstdio>
#include <iostream>
#include <CL/cl.hpp>
#include "clutil.h"

void CLkernelEwald(CLCTX cl);


static void setup_ewald(CLCTX ctx) {
    cl_int rc;
    ctx->queueEwald = clCreateCommandQueue(ctx->context->clContext, ctx->context->clDeviceId, 0, &rc);
    assert(rc == CL_SUCCESS);
    ctx->eventEwald = clCreateUserEvent(ctx->context->clContext, &rc);
    assert(rc == CL_SUCCESS);
    ctx->context->LxEwald = NULL;
    ctx->context->LyEwald = NULL;
    ctx->context->LxEwald = NULL;
    ctx->context->bhEwald = NULL;
    ctx->context->ewEwald = NULL;
    ctx->context->hxEwald = NULL;
    ctx->context->hyEwald = NULL;
    ctx->context->hzEwald = NULL;
    ctx->context->hCfac = NULL;
    ctx->context->hSfac = NULL;
    }


cl_program CL_compile(CLCTX cl, const char *src) {
    cl_int rc;
    cl_program program = clCreateProgramWithSource(cl->context->clContext, 1, &src, NULL, &rc);
    assert(rc == CL_SUCCESS);
    rc = clBuildProgram(program, 0, NULL,
	"-cl-strict-aliasing -cl-mad-enable -cl-denorms-are-zero -cl-fast-relaxed-math", NULL, NULL);

    if (rc != CL_SUCCESS) {
	size_t logSize;
	cl_build_status status;
	clGetProgramBuildInfo(program, cl->context->clDeviceId,
	    CL_PROGRAM_BUILD_STATUS, sizeof(cl_build_status), &status, NULL);
        clGetProgramBuildInfo(program, cl->context->clDeviceId, 
	    CL_PROGRAM_BUILD_LOG, 0, NULL, &logSize);
        char *programLog = (char*) calloc (logSize+1, sizeof(char));
        clGetProgramBuildInfo(program, cl->context->clDeviceId, 
                CL_PROGRAM_BUILD_LOG, logSize+1, programLog, NULL);
        printf("Build failed; error=%d, status=%d, programLog:nn%s", 
                rc, status, programLog);
        free(programLog);
	}
    assert(rc == CL_SUCCESS);

#if 0
    size_t bin_sz;
    rc = clGetProgramInfo(cl->context->programEwald, CL_PROGRAM_BINARY_SIZES, sizeof(size_t), &bin_sz, NULL);
 
    // Read binary (PTX file) to memory buffer
    unsigned char *bin = (unsigned char *)malloc(bin_sz);
    rc = clGetProgramInfo(cl->context->programEwald, CL_PROGRAM_BINARIES, sizeof(unsigned char *), &bin, NULL);
 
    // Save PTX to add_vectors_ocl.ptx
    FILE *fp = fopen("clewald.ptx", "wb");
    fwrite(bin, sizeof(char), bin_sz, fp);
    fclose(fp);
    free(bin);
 #endif

    return program;
    }


extern "C"
void *CL_create_context() {
    std::vector<cl::Platform> platforms;
    cl::Platform::get(&platforms);
    CLCONTEXT context = reinterpret_cast<CLCONTEXT>(malloc(sizeof(struct clContext)));
    assert(context != NULL);

    context->clContext = NULL;
    context->clDeviceId = NULL;

    if (platforms.size() > 0) {
	cl_int rc;
	rc = clGetDeviceIDs(platforms[0](),CL_DEVICE_TYPE_GPU,1,&context->clDeviceId,NULL);
	assert(rc == CL_SUCCESS);
	cl_context_properties properties[] = 
	    { CL_CONTEXT_PLATFORM, (cl_context_properties)(platforms[0])(), 0};
	context->clContext = clCreateContextFromType(properties, CL_DEVICE_TYPE_GPU, NULL, NULL, &rc);
	assert(rc == CL_SUCCESS);
	}
    return context;
    }


extern "C"
void *CL_initialize(void *vctx,int nCores,int iCore) {
    CLCONTEXT context = reinterpret_cast<CLCONTEXT>(vctx);
    cl_int rc;
    CLCTX ctx = reinterpret_cast<CLCTX>(malloc(sizeof(struct cl_ctx)));
    assert(ctx!=NULL);
    ctx->nCores = nCores;
    ctx->iCore = iCore;
    //ctx->epoch = 0;
    ctx->nWorkQueueSize = 0;
    ctx->nWorkQueueBusy = 0;
    ctx->wqCuda = ctx->wqFree = NULL;
    ctx->ewIn = NULL;
    ctx->ewt = NULL;

    ctx->context = context;
    CLkernelEwald(ctx);
    return ctx;
    }

extern "C"
void CL_SetQueueSize(void *vcl,int clSize, int inClBufSize, int outClBufSize) {
    CLCTX cl = reinterpret_cast<CLCTX>(vcl);
    cl_int rc;

    cl->inClBufSize = inClBufSize;
    cl->outClBufSize = outClBufSize;
    cl->nWorkQueueSize = clSize < 6 ? 6 : clSize;

    int hostBufSize = cl->inClBufSize > cl->outClBufSize ? cl->inClBufSize : cl->outClBufSize;
    while(clSize--) {
	CLwqNode *work = reinterpret_cast<CLwqNode *>(malloc(sizeof(CLwqNode)));

	work->clQueue = clCreateCommandQueue(cl->context->clContext, cl->context->clDeviceId, 0, &rc);
	assert(rc == CL_SUCCESS);
	work->clEvent = clCreateUserEvent(cl->context->clContext, &rc);
	assert(rc == CL_SUCCESS);

	// Pinned host memory buffers
	work->memOutCPU = clCreateBuffer(cl->context->clContext,CL_MEM_READ_ONLY|CL_MEM_ALLOC_HOST_PTR,hostBufSize,NULL,&rc);
	assert(rc == CL_SUCCESS);
	work->memInCPU  = clCreateBuffer(cl->context->clContext,CL_MEM_WRITE_ONLY|CL_MEM_ALLOC_HOST_PTR,hostBufSize,NULL,&rc);
	assert(rc == CL_SUCCESS);

	// GPU memory
	work->memInGPU  = clCreateBuffer(cl->context->clContext,CL_MEM_READ_ONLY,hostBufSize,NULL,&rc);
	assert(rc == CL_SUCCESS);
	work->memOutGPU = clCreateBuffer(cl->context->clContext,CL_MEM_WRITE_ONLY,hostBufSize,NULL,&rc);
	assert(rc == CL_SUCCESS);

	work->pHostBufToGPU = clEnqueueMapBuffer(work->clQueue,work->memOutCPU,CL_TRUE,CL_MAP_WRITE,0,hostBufSize,0,NULL,NULL,&rc);
	assert(rc == CL_SUCCESS);
	work->pHostBufFromGPU = clEnqueueMapBuffer(work->clQueue,work->memInCPU,CL_TRUE,CL_MAP_READ,0,hostBufSize,0,NULL,NULL,&rc);
	assert(rc == CL_SUCCESS);

        work->ctx = NULL;
        work->checkFcn = NULL;
        work->startTime = 0;
        work->next = cl->wqFree;
        cl->wqFree = work;
        }
    setup_ewald(cl);
    }

extern "C"
int CL_flushDone(void *vcl) {
    CLCTX cl = reinterpret_cast<CLCTX>(vcl);
    CLwqNode *work, **last = &cl->wqCuda;
    // Someone else has started the recovery process -- we just participate
//    if (cuda->epoch != recovery_epoch) {
//        CUDA_attempt_recovery(cuda,cudaSuccess);
//        }
    while( (work=*last) !=NULL ) {
	cl_int status;
	cl_int rc = clGetEventInfo(work->clEvent,CL_EVENT_COMMAND_EXECUTION_STATUS,sizeof(status),&status,NULL);
	if (rc!=CL_SUCCESS) std::clog << rc << " " << work << " " << work->clEvent << std::endl;
	assert(rc==CL_SUCCESS);
        if (status==CL_COMPLETE) {
            assert(work->checkFcn != NULL); /* Only one cudaSuccess per customer! */
            int rc = (*work->checkFcn)(work->ctx,work);
            assert(rc == 0);
            work->checkFcn = NULL; /* Make sure we don't do this multiple times! */
            *last = work->next;
            work->next = cl->wqFree;
            cl->wqFree = work;

            --cl->nWorkQueueBusy;
            continue;
            }
//        else if (rc!=cudaErrorNotReady) {
//            CUDA_attempt_recovery(cuda,rc);
//            }
//        else if (work->startTime != 0) {
//            double seconds = CUDA_getTime() - work->startTime;
//            if (seconds>=1.0) {
//                fprintf(stderr,"%s: cudaXxxxxQuery has returned cudaErrorNotReady for %f seconds\n",cuda->hostname,seconds);
//                work->startTime = 0;
//                CUDA_attempt_recovery(cuda,cudaErrorLaunchTimeout);
//                break;
//                }
//            }
        last = &work->next;
        }
    return cl->wqCuda != NULL;
    }



extern "C"
int CL_queue(void *vcl, void *ctx,
    int (*initWork)(void *vcl,void *ctx,void *work),
    int (*checkWork)(void *ctx,void *work),
    int (*doneWork)(void *ctx)) {
    CLCTX cl = reinterpret_cast<CLCTX>(vcl);
    CLwqNode *work;
    CL_flushDone(vcl);
    if (cl->wqFree == NULL || initWork==NULL) return 0;
    work = cl->wqFree;
    cl->wqFree = work->next;
    ++cl->nWorkQueueBusy;
    work->ctx = ctx;
    work->checkFcn = checkWork;
    work->initFcn = initWork;
    work->next = cl->wqCuda;
    cl->wqCuda = work;
//    work->startTime = CUDA_getTime();
    int rc = ((*work->initFcn)(vcl,work->ctx,work));
//    if ( rc != cudaSuccess) CUDA_attempt_recovery(cuda,rc);
    return 1;
    }
