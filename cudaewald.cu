#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <stdio.h>
#include "basetype.h"
#include "moments.h"
#include "cudautil.h"

#define cfOffset(i) ((i) + ((i)>>CUDA_LOG_NUM_BANKS))

#define MAX_TOTAL_REPLICAS (7*7*7)
#define SCAN_SIZE 512 // Must be larger than MAX_TOTAL_REPLICAS and a power of two.

__constant__ struct EwaldVariables ew;
__constant__ float hx[MAX_TOTAL_REPLICAS];
__constant__ float hy[MAX_TOTAL_REPLICAS];
__constant__ float hz[MAX_TOTAL_REPLICAS];
__constant__ float hCfac[MAX_TOTAL_REPLICAS];
__constant__ float hSfac[MAX_TOTAL_REPLICAS];


/*
** threadIdx.x: all work on the same particle -- this is the warp size, i.e., 32
** blockIdx.x:  different particles
**
** We are allowed 16 resident blocks so this corresponds to 512 threads per SM.
** This is fine because we are actually close to shared memory limited:
**   bValid[512] = 2052 bytes
**   We have 48k per SM so around 23 active thread blocks (but we use only 16)
**   compute 5.0: 32 thread blocks, but 64K shared gives us ~ 31 resident.
*/

__global__ void cudaEwald(double *inout) {
    __shared__ int bValid[cfOffset(SCAN_SIZE)]; // must be more than replicas
    double x, y, z;
    double r2,dir,dir2,a;
    double rx, ry, rz;
    double g0,g1,g2,g3,g4,g5,alphan;
    int nGrid= 2*ew.nEwReps + 1;
    int i000 = 3 + nGrid*(3 + 3*nGrid);
    int i, j, w, ix, iy, iz, bInHole;
    int pidx = blockIdx.x + gridDim.x*(blockIdx.y + gridDim.y*blockIdx.z);
    double tax = 0.0, tay = 0.0, taz = 0.0, tpot=0.0;
    int nP =  gridDim.x * gridDim.y * gridDim.z;

    double *X = inout;
    double *Y = inout + 1*nP;
    double *Z = inout + 2*nP;
    double *pPot = inout + 3*nP;

    rx = X[pidx] - ew.r[0];
    ry = Y[pidx] - ew.r[1];
    rz = Z[pidx] - ew.r[2];

    // Step 1: determine if we need to process this element
    // There are 7 x 7 x 7 replicas, and we need to process only ~ 25% of them.
    // We specifically let the CPU handle 0,0,0 because it can be special
    for( i=threadIdx.x; i<MAX_TOTAL_REPLICAS; i+= blockDim.x) {
	// Calculate the indexes for the replicas from the thread index
	// The range is [-nEwReps,nEwReps], normally [3,3]
	iz = i / (nGrid*nGrid);
	j = i - iz*(nGrid*nGrid);
	iy = j / nGrid;
	ix = j - iy*nGrid;
	ix -= ew.nEwReps;
	iy -= ew.nEwReps;
	iz -= ew.nEwReps;
	// From that we can calculate the distance
	x = rx + ew.Lbox * ix;
	y = ry + ew.Lbox * iy;
	z = rz + ew.Lbox * iz;
	r2 = x*x + y*y + z*z;
	bInHole = (abs(ix) <= ew.nReps && abs(iy) <= ew.nReps && abs(iz) <= ew.nReps);
	// We omit 0,0,0 because it can be special. We let the CPU handle it.
	bValid[1+cfOffset(i)] = i!=i000 && (r2 < ew.fEwCut2 || bInHole);
	}
    for(; i<SCAN_SIZE; i+= blockDim.x) bValid[1+cfOffset(i)] = 0;

    // Step 2: scan to get the valid elements
    for(w=2; w<=SCAN_SIZE; w*=2) { // scan up
	for( i=threadIdx.x*w; i<SCAN_SIZE; i+= w*blockDim.x ) {
	    bValid[1+cfOffset(i+w-1)] += bValid[1+cfOffset(i+w/2-1)];
	    }
	}
    if (threadIdx.x==0) bValid[1+cfOffset(SCAN_SIZE-1)] = 0;
    for(w=SCAN_SIZE; w>1; w/=2) { // scan down
	for( i=threadIdx.x*w; i<SCAN_SIZE; i+= w*blockDim.x ) {
            j = bValid[1+cfOffset(i+w/2-1)];
            bValid[1+cfOffset(i+w/2-1)] = bValid[1+cfOffset(i+w-1)];
            bValid[1+cfOffset(i+w-1)] += j;
	    }
	}

    int nCompute = bValid[1+cfOffset(SCAN_SIZE-1)];

    // Flatten the array - avoids bank conflicts in the next step
    for( i=threadIdx.x; i<SCAN_SIZE; i+= blockDim.x) {
	bValid[1+i] = bValid[1+cfOffset(i)];
	}

    // Step 3: create a list of work
    if (threadIdx.x==0) bValid[0] = 0;
    for( i=threadIdx.x; i<SCAN_SIZE; i+= blockDim.x) {
	if ( bValid[1+i-1] != bValid[1+i]) {
	    bValid[bValid[1+i]-1] = i-1;
	    }
	}

    // At this point we have around 80 some elements (with the default parameters)
    // This means we will loop 3 times. We have cut the work by more than a factor of four.
    for( i=threadIdx.x; i<nCompute; i+= blockDim.x) {
	j = bValid[i];

	iz = j / (nGrid*nGrid);
	j = j - iz*(nGrid*nGrid);
	iy = j / nGrid;
	ix = j - iy*nGrid;
	ix -= ew.nEwReps;
	iy -= ew.nEwReps;
	iz -= ew.nEwReps;

	bInHole = (abs(ix) <= ew.nReps && abs(iy) <= ew.nReps && abs(iz) <= ew.nReps);

	x = rx + ew.Lbox * ix;
	y = ry + ew.Lbox * iy;
	z = rz + ew.Lbox * iz;
	r2 = x*x + y*y + z*z;

	//dir = 1/sqrt(r2);
	dir = rsqrt(r2);
	dir2 = dir*dir;
	a = exp(-r2*ew.alpha2);
	a *= ew.ka*dir2;
	if (bInHole) g0 = -erf(ew.alpha*r2*dir);
	else 	 g0 = erfc(ew.alpha*r2*dir);
	g0 *= dir;
	g1 = g0*dir2 + a;
	alphan = 2*ew.alpha2;
	g2 = 3*g1*dir2 + alphan*a;
	alphan *= 2*ew.alpha2;
	g3 = 5*g2*dir2 + alphan*a;
	alphan *= 2*ew.alpha2;
	g4 = 7*g3*dir2 + alphan*a;
	alphan *= 2*ew.alpha2;
	g5 = 9*g4*dir2 + alphan*a;

	double onethird = 1.0/3.0;
	double xx,xxx,xxy,xxz,yy,yyy,yyz,xyy,zz,zzz,xzz,yzz,xy,xyz,xz,yz;
	double Qta,Q4mirx,Q4miry,Q4mirz,Q4mir,Q4x,Q4y,Q4z;
	double Q3mirx,Q3miry,Q3mirz,Q3mir,Q2mirx,Q2miry,Q2mirz,Q2mir;

	xx = 0.5*x*x;
	xxx = onethird*xx*x;
	xxy = xx*y;
	xxz = xx*z;
	yy = 0.5*y*y;
	yyy = onethird*yy*y;
	xyy = yy*x;
	yyz = yy*z;
	zz = 0.5*z*z;
	zzz = onethird*zz*z;
	xzz = zz*x;
	yzz = zz*y;
	xy = x*y;
	xyz = xy*z;
	xz = x*z;
	yz = y*z;
	Q2mirx = ew.mom.xx*x + ew.mom.xy*y + ew.mom.xz*z;
	Q2miry = ew.mom.xy*x + ew.mom.yy*y + ew.mom.yz*z;
	Q2mirz = ew.mom.xz*x + ew.mom.yz*y + ew.mom.zz*z;
	Q3mirx = ew.mom.xxx*xx + ew.mom.xxy*xy + ew.mom.xxz*xz + ew.mom.xyy*yy + ew.mom.xyz*yz + ew.mom.xzz*zz;
	Q3miry = ew.mom.xxy*xx + ew.mom.xyy*xy + ew.mom.xyz*xz + ew.mom.yyy*yy + ew.mom.yyz*yz + ew.mom.yzz*zz;
	Q3mirz = ew.mom.xxz*xx + ew.mom.xyz*xy + ew.mom.xzz*xz + ew.mom.yyz*yy + ew.mom.yzz*yz + ew.mom.zzz*zz;
	Q4mirx = ew.mom.xxxx*xxx + ew.mom.xxxy*xxy + ew.mom.xxxz*xxz + ew.mom.xxyy*xyy + ew.mom.xxyz*xyz +
	    ew.mom.xxzz*xzz + ew.mom.xyyy*yyy + ew.mom.xyyz*yyz + ew.mom.xyzz*yzz + ew.mom.xzzz*zzz;
	Q4miry = ew.mom.xxxy*xxx + ew.mom.xxyy*xxy + ew.mom.xxyz*xxz + ew.mom.xyyy*xyy + ew.mom.xyyz*xyz +
	    ew.mom.xyzz*xzz + ew.mom.yyyy*yyy + ew.mom.yyyz*yyz + ew.mom.yyzz*yzz + ew.mom.yzzz*zzz;
	Q4mirz = ew.mom.xxxz*xxx + ew.mom.xxyz*xxy + ew.mom.xxzz*xxz + ew.mom.xyyz*xyy + ew.mom.xyzz*xyz +
	    ew.mom.xzzz*xzz + ew.mom.yyyz*yyy + ew.mom.yyzz*yyz + ew.mom.yzzz*yzz + ew.mom.zzzz*zzz;
	Q4x = ew.Q4xx*x + ew.Q4xy*y + ew.Q4xz*z;
	Q4y = ew.Q4xy*x + ew.Q4yy*y + ew.Q4yz*z;
	Q4z = ew.Q4xz*x + ew.Q4yz*y + ew.Q4zz*z;
	Q2mir = 0.5*(Q2mirx*x + Q2miry*y + Q2mirz*z) - (ew.Q3x*x + ew.Q3y*y + ew.Q3z*z) + ew.Q4;
	Q3mir = onethird*(Q3mirx*x + Q3miry*y + Q3mirz*z) - 0.5*(Q4x*x + Q4y*y + Q4z*z);
	Q4mir = 0.25*(Q4mirx*x + Q4miry*y + Q4mirz*z);
	Qta = g1*ew.mom.m - g2*ew.Q2 + g3*Q2mir + g4*Q3mir + g5*Q4mir;
	tpot -= g0*ew.mom.m - g1*ew.Q2 + g2*Q2mir + g3*Q3mir + g4*Q4mir;
	tax += g2*(Q2mirx - ew.Q3x) + g3*(Q3mirx - Q4x) + g4*Q4mirx - x*Qta;
	tay += g2*(Q2miry - ew.Q3y) + g3*(Q3miry - Q4y) + g4*Q4miry - y*Qta;
	taz += g2*(Q2mirz - ew.Q3z) + g3*(Q3mirz - Q4z) + g4*Q4mirz - z*Qta;
	}

    // the H-Loop
    float fx=rx, fy=ry, fz=rz;
    for( i=threadIdx.x; i<ew.nEwhLoop; i+= blockDim.x) {
	float hdotx,s,c,t;
	hdotx = hx[i]*fx + hy[i]*fy + hz[i]*fz;
	sincosf(hdotx,&s,&c);
	tpot += hCfac[i]*c + hSfac[i]*s;
	t = hCfac[i]*s - hSfac[i]*c;
	tax += hx[i]*t;
	tay += hy[i]*t;
	taz += hz[i]*t;
	}

#define R1(S,T,O) (S)[threadIdx.x] = T = T + (S)[threadIdx.x + O]
#define R(S,T) R1(S,T,16); R1(S,T,8); R1(S,T,4); R1(S,T,2); R1(S,T,1);

    volatile float * v = (float *)bValid; // We can use this buffer now
    v[threadIdx.x] = tax;
    v[threadIdx.x+32] = tay;
    v[threadIdx.x+64] = taz;
    v[threadIdx.x+96] = tpot;
    if (threadIdx.x < 16) {
	R(v,tax); R(v+32,tay); R(v+64,taz); R(v+96,tpot);
	if (threadIdx.x==0) {
	    pPot[pidx] = tpot;
	    X[pidx] = tax;
	    Y[pidx] = tay;
	    Z[pidx] = taz;
	    }
	}
    }

extern "C"
void cudaEwaldInit(struct EwaldVariables *ewIn, EwaldTable *ewt ) {
    CUDA_CHECK(cudaMemcpyToSymbol,(ew,ewIn,sizeof(ew),0,cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpyToSymbol,(hx,ewt->hx.f,sizeof(float)*ewIn->nEwhLoop,0,cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpyToSymbol,(hy,ewt->hy.f,sizeof(float)*ewIn->nEwhLoop,0,cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpyToSymbol,(hz,ewt->hz.f,sizeof(float)*ewIn->nEwhLoop,0,cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpyToSymbol,(hCfac,ewt->hCfac.f,sizeof(float)*ewIn->nEwhLoop,0,cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpyToSymbol,(hSfac,ewt->hSfac.f,sizeof(float)*ewIn->nEwhLoop,0,cudaMemcpyHostToDevice));
    }

extern "C"
int CUDAinitWorkEwald( void *ve, void *vwork ) {
    workEwald *e = reinterpret_cast<workEwald *>(ve);
    CUDAwqNode *work = reinterpret_cast<CUDAwqNode *>(vwork);
    //CUDACTX ctx = reinterpret_cast<CUDACTX>(e->cudaCtx);
    double *pHostBuf = reinterpret_cast<double *>(work->pHostBuf);
    double *pCudaBuf = reinterpret_cast<double *>(work->pCudaBuf);
    PARTICLE *p;
    int nx, ny, nz, i;

    // cuda global arrays
    nx = e->nP < 65536 ? e->nP : 65536;
    ny = e->nP / nx;
    if (ny==0) ny = 1;
    else if (ny>=65536) ny=65536;
    nz = e->nP / (nx*ny);
    if (nz==0) nz = 1;
    assert(nx*ny*nz == e->nP); /* We do powers of two only */

    dim3 dimBlock( 32, 1 );
    dim3 dimGrid( nx, ny,nz );

    for(i=0; i<e->nP; ++i) {
	p = e->pPart[i];
	pHostBuf[i + 0*e->nP] = p->r[0];
	pHostBuf[i + 1*e->nP] = p->r[1];
	pHostBuf[i + 2*e->nP] = p->r[2];
	}

    // copy data directly to device memory
    CUDA_CHECK(cudaMemcpyAsync,(pCudaBuf, pHostBuf, e->nP*3*sizeof(double),
	    cudaMemcpyHostToDevice, work->stream));
//    CUDA_CHECK(cudaMemcpy,(pCudaBuf, pHostBuf, e->nP*3*sizeof(double),
//	    cudaMemcpyHostToDevice));
    cudaEwald<<<dimGrid, dimBlock, 0, work->stream>>>(pCudaBuf);
//    cudaEwald<<<dimGrid, dimBlock>>>(pCudaBuf);
    CUDA_CHECK(cudaMemcpyAsync,(pHostBuf, pCudaBuf, e->nP*4*sizeof(double),
            cudaMemcpyDeviceToHost, work->stream));
//    CUDA_CHECK(cudaMemcpy,(pHostBuf, pCudaBuf, e->nP*4*sizeof(double),
//            cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaEventRecord,(work->event,work->stream));

    return 1;
    }

extern "C"
void pkdAccumulateCUDA(void * pkd,int nP,PARTICLE **pPart,double *pax,double *pay,double *paz,double *pot);


extern "C"
int CUDAcheckWorkEwald( void *ve, void *vwork ) {
    workEwald *e = reinterpret_cast<workEwald *>(ve);
    CUDAwqNode *work = reinterpret_cast<CUDAwqNode *>(vwork);
    double *pHostBuf = reinterpret_cast<double *>(work->pHostBuf);

    pkdAccumulateCUDA(e->pkd,e->nP,e->pPart,
	pHostBuf + 0*e->nP, pHostBuf + 1*e->nP,
	pHostBuf + 2*e->nP, pHostBuf + 3*e->nP);
    free(e->pPart);
    free(e);
    return 0;

    }
