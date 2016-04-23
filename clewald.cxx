#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "clutil.h"
#include <iostream>

#define ALIGN 64
#define MASK (ALIGN-1)

const char *ewald_kernel_src =
    "typedef struct momComplete {\n"
    "    double m;\n"
    "    double xx,yy,xy,xz,yz;\n"
    "    double xxx,xyy,xxy,yyy,xxz,yyz,xyz;\n"
    "    double xxxx,xyyy,xxxy,yyyy,xxxz,yyyz,xxyy,xxyz,xyyz;\n"
    "    double zz;\n"
    "    double xzz,yzz,zzz;\n"
    "    double xxzz,xyzz,xzzz,yyzz,yzzz,zzzz;\n"
    "    } MOMC;\n"
    "typedef struct {\n"
    "    double r[3];\n"
    "    MOMC mom;\n"
    "    double fEwCut2,fInner2,alpha,ialpha,alpha2,k1,ka,Lbox;\n"
    "    double Q4xx,Q4xy,Q4xz,Q4yy,Q4yz,Q4zz,Q4,Q3x,Q3y,Q3z,Q2;\n"
    "    int nMaxEwhLoop;\n"
    "    int nEwLoopInner, nEwhLoop;\n"
    "    int nReps,nEwReps;\n"
    "    } EwaldVariables;\n"
    "__kernel void ewald(__constant EwaldVariables * ew,\n"
    "    __constant float *hx,\n"
    "    __constant float *hy,\n"
    "    __constant float *hz,\n"
    "    __constant float *hCfac,\n"
    "    __constant float *hSfac,\n"
    "    int nAlign,\n"
    "    __global double *dIn,\n"
    "    __global double *dOut)\n"
    "    {\n"
    "    __global double *X = dIn + 0*nAlign;\n"
    "    __global double *Y = dIn + 1*nAlign;\n"
    "    __global double *Z = dIn + 2*nAlign;\n"
    "    __global double *Xout = dOut + 0*nAlign;\n"
    "    __global double *Yout = dOut + 1*nAlign;\n"
    "    __global double *Zout = dOut + 2*nAlign;\n"
    "    __global double *pPot = dOut + 3*nAlign;\n"
    "    __global double *pdFlop = dOut + 4*nAlign;\n"
    "    double g0,g1,g2,g3,g4,g5;\n"
    "    int i, ix, iy, iz, bInHole;\n"
    "    int pidx = get_global_id(0);\n"
    "    double tax = 0.0, tay = 0.0, taz = 0.0, tpot=0.0, dFlop=0.0;\n"
    "    const double rx = X[pidx] - ew->r[0];\n"
    "    const double ry = Y[pidx] - ew->r[1];\n"
    "    const double rz = Z[pidx] - ew->r[2];\n"
    "    for(ix=-3; ix<=3; ++ix) {\n"
    "        for(iy=-3; iy<=3; ++iy) {\n"
    "            for(iz=-3; iz<=3; ++iz) {\n"
    "                bInHole = (abs(ix) <= ew->nReps && abs(iy) <= ew->nReps && abs(iz) <= ew->nReps);\n"
    "                const double x = rx + ew->Lbox * ix;\n"
    "                const double y = ry + ew->Lbox * iy;\n"
    "                const double z = rz + ew->Lbox * iz;\n"
    "                double r2 = x*x + y*y + z*z;\n"
    "                if (r2 >= ew->fEwCut2 && !bInHole) continue;\n"
    "                if (r2 < ew->fInner2) {\n"
    "                    double alphan = ew->ka;\n"
    "                    r2 *= ew->alpha2;\n"
    "                    g0 = alphan*((1.0/3.0)*r2 - 1.0);\n"
    "                    alphan *= 2*ew->alpha2;\n"
    "                    g1 = alphan*((1.0/5.0)*r2 - (1.0/3.0));\n"
    "                    alphan *= 2*ew->alpha2;\n"
    "                    g2 = alphan*((1.0/7.0)*r2 - (1.0/5.0));\n"
    "                    alphan *= 2*ew->alpha2;\n"
    "                    g3 = alphan*((1.0/9.0)*r2 - (1.0/7.0));\n"
    "                    alphan *= 2*ew->alpha2;\n"
    "                    g4 = alphan*((1.0/11.0)*r2 - (1.0/9.0));\n"
    "                    alphan *= 2*ew->alpha2;\n"
    "                    g5 = alphan*((1.0/13.0)*r2 - (1.0/11.0));\n"
    "                    }\n"
    "                else {\n"
    "                    double dir = rsqrt(r2);\n"
    "                    double dir2 = dir*dir;\n"
    "                    double a = exp(-r2*ew->alpha2);\n"
    "                    a *= ew->ka*dir2;\n"
    "                    if (bInHole) g0 = -erf(ew->alpha*r2*dir);\n"
    "                    else         g0 = erfc(ew->alpha*r2*dir);\n"
    "                    g0 *= dir;\n"
    "                    g1 = g0*dir2 + a;\n"
    "                    double alphan = 2*ew->alpha2;\n"
    "                    g2 = 3*g1*dir2 + alphan*a;\n"
    "                    alphan *= 2*ew->alpha2;\n"
    "                    g3 = 5*g2*dir2 + alphan*a;\n"
    "                    alphan *= 2*ew->alpha2;\n"
    "                    g4 = 7*g3*dir2 + alphan*a;\n"
    "                    alphan *= 2*ew->alpha2;\n"
    "                    g5 = 9*g4*dir2 + alphan*a;\n"
    "                    }\n"
    "                double onethird = 1.0/3.0;\n"
    "                double Qta,Q4mirx,Q4miry,Q4mirz,Q4mir,Q4x,Q4y,Q4z;\n"
    "                double Q3mirx,Q3miry,Q3mirz,Q3mir,Q2mirx,Q2miry,Q2mirz,Q2mir;\n"
    "                tpot -= g0*ew->mom.m - g1*ew->Q2;"
    "                Qta = g1*ew->mom.m - g2*ew->Q2;"

    "                const double xx = 0.5*x*x;\n"
    "                const double xxx = onethird*xx*x;\n"
    "                const double xxy = xx*y;\n"
    "                const double xxz = xx*z;\n"
    "                const double yy = 0.5*y*y;\n"
    "                const double yyy = onethird*yy*y;\n"
    "                const double xyy = yy*x;\n"
    "                const double yyz = yy*z;\n"
    "                const double zz = 0.5*z*z;\n"
    "                const double zzz = onethird*zz*z;\n"
    "                const double xzz = zz*x;\n"
    "                const double yzz = zz*y;\n"
    "                const double xy = x*y;\n"
    "                const double xyz = xy*z;\n"
    "                const double xz = x*z;\n"
    "                const double yz = y*z;\n"

    "                Q2mir = -(ew->Q3x*x + ew->Q3y*y + ew->Q3z*z) + ew->Q4;\n"
    "                Q2mirx = ew->mom.xx*x + ew->mom.xy*y + ew->mom.xz*z;\n"
    "                tax += g2*(Q2mirx - ew->Q3x);"
    "                Q2miry = ew->mom.xy*x + ew->mom.yy*y + ew->mom.yz*z;\n"
    "                tay += g2*(Q2miry - ew->Q3y);"
    "                Q2mirz = ew->mom.xz*x + ew->mom.yz*y + ew->mom.zz*z;\n"
    "                taz += g2*(Q2mirz - ew->Q3z);"
    "                Q2mir += 0.5*(0.5*Q2mirx*x + 0.5*Q2miry*y + Q2mirz*z);"
    "                Qta += g3*Q2mir; tpot -= g2*Q2mir;"

    "                Q4x = ew->Q4xx*x + ew->Q4xy*y + ew->Q4xz*z;\n"
    "                Q4y = ew->Q4xy*x + ew->Q4yy*y + ew->Q4yz*z;\n"
    "                Q4z = ew->Q4xz*x + ew->Q4yz*y + ew->Q4zz*z;\n"
    "                Q3mir = -0.5*(Q4x*x + Q4y*y + Q4z*z);\n"

    "                Q3mirx = ew->mom.xxx*xx + ew->mom.xxy*xy + ew->mom.xxz*xz + ew->mom.xyy*yy + ew->mom.xyz*yz + ew->mom.xzz*zz;\n"
    "                tax += g3*(Q3mirx - Q4x);"
    "                Q3miry = ew->mom.xxy*xx + ew->mom.xyy*xy + ew->mom.xyz*xz + ew->mom.yyy*yy + ew->mom.yyz*yz + ew->mom.yzz*zz;\n"
    "                tay += g3*(Q3miry - Q4y);"
    "                Q3mirz = ew->mom.xxz*xx + ew->mom.xyz*xy + ew->mom.xzz*xz + ew->mom.yyz*yy + ew->mom.yzz*yz + ew->mom.zzz*zz;\n"
    "                taz += g3*(Q3mirz - Q4z);"
    "                Q3mir += onethird*(Q3mirx*x+Q3miry*y+Q3mirz*z);"
    "                Qta += g4*Q3mir; tpot -= g3*Q3mir;"

    "                Q4mirx = ew->mom.xxxx*xxx + ew->mom.xxxy*xxy + ew->mom.xxxz*xxz + ew->mom.xxyy*xyy + ew->mom.xxyz*xyz +\n"
    "                    ew->mom.xxzz*xzz + ew->mom.xyyy*yyy + ew->mom.xyyz*yyz + ew->mom.xyzz*yzz + ew->mom.xzzz*zzz;\n"
    "                tax += g4*Q4mirx;"
    "                Q4miry = ew->mom.xxxy*xxx + ew->mom.xxyy*xxy + ew->mom.xxyz*xxz + ew->mom.xyyy*xyy + ew->mom.xyyz*xyz +\n"
    "                    ew->mom.xyzz*xzz + ew->mom.yyyy*yyy + ew->mom.yyyz*yyz + ew->mom.yyzz*yzz + ew->mom.yzzz*zzz;\n"
    "                tay += g4*Q4miry;"
    "                Q4mirz = ew->mom.xxxz*xxx + ew->mom.xxyz*xxy + ew->mom.xxzz*xxz + ew->mom.xyyz*xyy + ew->mom.xyzz*xyz +\n"
    "                    ew->mom.xzzz*xzz + ew->mom.yyyz*yyy + ew->mom.yyzz*yyz + ew->mom.yzzz*yzz + ew->mom.zzzz*zzz;\n"
    "                taz += g4*Q4mirz;"
    "                Q4mir = 0.25*(Q4mirx*x+Q4miry*y+Q4mirz*z);"

    "                Qta += g5*Q4mir; tpot -= g4*Q4mir;"
    "                tax -= x*Qta;\n"
    "                tay -= y*Qta;\n"
    "                taz -= z*Qta;\n"
    "                dFlop += 386;\n"
    "                }\n"
    "            }\n"
    "        }\n"
    "    float fx=rx, fy=ry, fz=rz;\n"
    "    float fax=0, fay=0, faz=0, fpot=0;\n"
    "    for( i=0; i<ew->nEwhLoop; ++i) {\n"
    "        float hdotx,s,c,t;\n"
    "	     hdotx = hx[i]*fx + hy[i]*fy + hz[i]*fz;\n"
    "        s = sincos(hdotx,&c);\n"
    "        fpot += hCfac[i]*c + hSfac[i]*s;\n"
    "        t = hCfac[i]*s - hSfac[i]*c;\n"
    "        fax += hx[i]*t;\n"
    "        fay += hy[i]*t;\n"
    "        faz += hz[i]*t;\n"
    "        }\n"
    "    Xout[pidx] = tax + fax;\n"
    "    Yout[pidx] = tay + fay;\n"
    "    Zout[pidx] = taz + faz;\n"
    "    pPot[pidx] = tpot + fpot;\n"
    "    pdFlop[pidx] = dFlop;\n"
    "    return;\n"
    "    }\n";


#include <stdio.h>
/* Compile/create ewald kernel */
void CLkernelEwald(CLCTX cl) {
    cl_int rc;
    cl->context->programEwald = clCreateProgramWithSource(cl->context->clContext, 1, &ewald_kernel_src, NULL, &rc);
    assert(rc == CL_SUCCESS);
    rc = clBuildProgram(cl->context->programEwald, 0, NULL,
	"-cl-strict-aliasing -cl-mad-enable -cl-denorms-are-zero -cl-fast-relaxed-math", NULL, NULL);
    assert(rc == CL_SUCCESS);
    cl->kernelEwald = clCreateKernel(cl->context->programEwald, "ewald", &rc);
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
    }

cl_int cl_setup_ewald(CLCTX cl) {
    if (cl->ewIn && cl->ewt) {
        cl_int rc, status;
//        double start = CUDA_getTime();

	if (cl->context->hxEwald==NULL) {
	    cl->context->ewEwald = clCreateBuffer(cl->context->clContext,CL_MEM_READ_ONLY,sizeof(struct EwaldVariables),NULL,&rc);
	    assert(rc == CL_SUCCESS);
	    cl->context->hxEwald = clCreateBuffer(cl->context->clContext,CL_MEM_READ_ONLY,sizeof(float)*cl->ewIn->nEwhLoop,NULL,&rc);
	    assert(rc == CL_SUCCESS);
	    cl->context->hyEwald = clCreateBuffer(cl->context->clContext,CL_MEM_READ_ONLY,sizeof(float)*cl->ewIn->nEwhLoop,NULL,&rc);
	    assert(rc == CL_SUCCESS);
	    cl->context->hzEwald = clCreateBuffer(cl->context->clContext,CL_MEM_READ_ONLY,sizeof(float)*cl->ewIn->nEwhLoop,NULL,&rc);
	    assert(rc == CL_SUCCESS);
	    cl->context->hCfac = clCreateBuffer(cl->context->clContext,CL_MEM_READ_ONLY,sizeof(float)*cl->ewIn->nEwhLoop,NULL,&rc);
	    assert(rc == CL_SUCCESS);
	    cl->context->hSfac = clCreateBuffer(cl->context->clContext,CL_MEM_READ_ONLY,sizeof(float)*cl->ewIn->nEwhLoop,NULL,&rc);
	    assert(rc == CL_SUCCESS);
	    }

	rc = clEnqueueWriteBuffer(cl->queueEwald,cl->context->ewEwald,CL_FALSE,0,sizeof(struct EwaldVariables),cl->ewIn,0,NULL,NULL);
	assert(rc == CL_SUCCESS);
	rc = clEnqueueWriteBuffer(cl->queueEwald,cl->context->hxEwald,CL_FALSE,0,sizeof(float)*cl->ewIn->nEwhLoop,cl->ewt->hx.f,0,NULL,NULL);
	assert(rc == CL_SUCCESS);
	rc = clEnqueueWriteBuffer(cl->queueEwald,cl->context->hyEwald,CL_FALSE,0,sizeof(float)*cl->ewIn->nEwhLoop,cl->ewt->hy.f,0,NULL,NULL);
	assert(rc == CL_SUCCESS);
	rc = clEnqueueWriteBuffer(cl->queueEwald,cl->context->hzEwald,CL_FALSE,0,sizeof(float)*cl->ewIn->nEwhLoop,cl->ewt->hz.f,0,NULL,NULL);
	assert(rc == CL_SUCCESS);
	rc = clEnqueueWriteBuffer(cl->queueEwald,cl->context->hCfac,CL_FALSE,0,sizeof(float)*cl->ewIn->nEwhLoop,cl->ewt->hCfac.f,0,NULL,NULL);
	assert(rc == CL_SUCCESS);
	rc = clEnqueueWriteBuffer(cl->queueEwald,cl->context->hSfac,CL_FALSE,0,sizeof(float)*cl->ewIn->nEwhLoop,cl->ewt->hSfac.f,0,NULL,&cl->eventEwald);
	assert(rc == CL_SUCCESS);
        do {
	    rc = clGetEventInfo(cl->eventEwald,CL_EVENT_COMMAND_EXECUTION_STATUS,sizeof(status),&status,NULL);
	    if (rc!=CL_SUCCESS) std::cout << rc << std::endl;
	    assert(rc==CL_SUCCESS);
//            if (CUDA_getTime() - start > 1.0) {
//                return cudaErrorLaunchTimeout;
//                }
            } while (status!=CL_COMPLETE);
	return CL_SUCCESS;
        }
    return CL_SUCCESS;
    }


extern "C"
void pkdAccumulateCUDA(void * pkd,workEwald *we,double *pax,double *pay,double *paz,double *pot,double *pdFlop);

extern "C"
int CLcheckWorkEwald( void *ve, void *vwork ) {
    workEwald *e = reinterpret_cast<workEwald *>(ve);
    CLwqNode *work = reinterpret_cast<CLwqNode *>(vwork);
    double *pHostBuf = reinterpret_cast<double *>(work->pHostBufFromGPU);
    double *X, *Y, *Z, *pPot, *pdFlop;
    int align;
    cl_int rc;


    align = (e->nP+MASK)&~MASK; /* As above! Warp align the memory buffers */

//    double *pHostBuf = reinterpret_cast<double *>(clEnqueueMapBuffer(work->clQueue,work->memInCPU,CL_TRUE,CL_MAP_READ,0,align*5*sizeof(double),0,NULL,NULL,&rc));
//    assert(rc == CL_SUCCESS);



    X       = pHostBuf + 0*align;
    Y       = pHostBuf + 1*align;
    Z       = pHostBuf + 2*align;
    pPot    = pHostBuf + 3*align;
    pdFlop  = pHostBuf + 4*align;
    pkdAccumulateCUDA(e->pkd,e,X,Y,Z,pPot,pdFlop);

//    clEnqueueUnmapMemObject(work->clQueue,work->memInCPU,pHostBuf,0,NULL,NULL);


    free(e->ppWorkPart);
    free(e->piWorkPart);
    free(e);
    return 0;

    }


extern "C"
void clEwaldInit(void *clCtx, struct EwaldVariables *ewIn, EwaldTable *ewt ) {
    CLCTX cl = reinterpret_cast<CLCTX>(clCtx);
    cl->ewIn = ewIn;
    cl->ewt = ewt;
    if (cl->iCore==0) {
	cl_setup_ewald(cl);
	}
    }

extern "C"
int CLinitWorkEwald( void *vcl, void *ve, void *vwork ) {
    CLCTX cl = reinterpret_cast<CLCTX>(vcl);
    workEwald *e = reinterpret_cast<workEwald *>(ve);
    CLwqNode *work = reinterpret_cast<CLwqNode *>(vwork);
    double *pHostBufToGPU    = reinterpret_cast<double *>(work->pHostBufToGPU);
    double *X, *Y, *Z;
    int align, i;
    cl_int rc;

    align = (e->nP+MASK)&~MASK; /* Warp align the memory buffers */

    X       = pHostBufToGPU + 0*align;
    Y       = pHostBufToGPU + 1*align;
    Z       = pHostBufToGPU + 2*align;
    for(i=0; i<e->nP; ++i) {
        const workParticle *wp = e->ppWorkPart[i];
	const int wi = e->piWorkPart[i];
	const PINFOIN *in = &wp->pInfoIn[wi];
	X[i] = wp->c[0] + in->r[0];
	Y[i] = wp->c[1] + in->r[1];
	Z[i] = wp->c[2] + in->r[2];
	}
    for(i;i<align;++i) X[i]=Y[i]=Z[i] = 100;
    rc = clEnqueueWriteBuffer(work->clQueue,work->memInGPU,CL_FALSE,0,align*3*sizeof(double),work->pHostBufToGPU,0,NULL,NULL);
    assert(rc == CL_SUCCESS);

    clSetKernelArg(cl->kernelEwald, 0, sizeof(cl_mem), (void *)&cl->context->ewEwald);
    clSetKernelArg(cl->kernelEwald, 1, sizeof(cl_mem), (void *)&cl->context->hxEwald);
    clSetKernelArg(cl->kernelEwald, 2, sizeof(cl_mem), (void *)&cl->context->hyEwald);
    clSetKernelArg(cl->kernelEwald, 3, sizeof(cl_mem), (void *)&cl->context->hzEwald);
    clSetKernelArg(cl->kernelEwald, 4, sizeof(cl_mem), (void *)&cl->context->hCfac);
    clSetKernelArg(cl->kernelEwald, 5, sizeof(cl_mem), (void *)&cl->context->hSfac);
    clSetKernelArg(cl->kernelEwald, 6, sizeof(int),    (void *)&align);
    clSetKernelArg(cl->kernelEwald, 7, sizeof(cl_mem), (void *)&work->memInGPU);
    clSetKernelArg(cl->kernelEwald, 8, sizeof(cl_mem), (void *)&work->memOutGPU);

    size_t global_work_size[2] = {align};
    size_t local_work_size[2] = {ALIGN};
    rc = clEnqueueNDRangeKernel(work->clQueue, cl->kernelEwald, 1, NULL, global_work_size, local_work_size, 0, NULL, NULL);
    assert(rc == CL_SUCCESS);

    rc = clEnqueueReadBuffer(work->clQueue,work->memOutGPU,CL_FALSE,0,align*5*sizeof(double),work->pHostBufFromGPU,0,NULL,&work->clEvent);
    assert(rc == CL_SUCCESS);

    return CL_SUCCESS;
    }
