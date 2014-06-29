/* -*- Mode: C++; tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <stdio.h>
#include "cudautil.h"

#define SYNC_RATE 16  // Must be: 1, 2, 4, 8, 16
#define WIDTH 32

#define TB_THREADS 128
#define WARPS 4

#define PP_WU 128
#define PC_WU 32

#include "ilp.h"
#include "ilc.h"
#include "basetype.h"

/*
** The following are the basically the same as ILP_BLK and ILC_BLK,
** but we need to be able to alter their sizes.
*/

template<int n>
struct ilpBlk {
    float dx[n], dy[n], dz[n];    /* Offset from ilp->cx, cy, cz */
    float m[n];             /* Mass */
    float fourh2[n];        /* Softening: calculated */
    };

template<int n>
struct ilcBlk {
    float dx[n],dy[n],dz[n];
    float xxxx[n],xxxy[n],xxxz[n],xxyz[n],xxyy[n],yyyz[n],xyyz[n],xyyy[n],yyyy[n];
    float xxx[n],xyy[n],xxy[n],yyy[n],xxz[n],yyz[n],xyz[n];
    float xx[n],xy[n],xz[n],yy[n],yz[n];
    float x[n],y[n],z[n];
    float m[n],u[n];
    };

// One of these entries for each interaction block
struct ppWorkUnit {
    int nP;   // Number of particles
    int iP;   // Index of first particle
    int nI;   // Number of interactions in the block
    int iO;   // Index of the output block
    };


struct ppInput {
    float dx, dy, dz;
    float ax, ay, az;
    float fSoft2;
    float dImaga;
    };

/* Each thread block outputs this for each particle */
struct __align__(32) ppResult {
    float ax;
    float ay;
    float az;
    float fPot;
    float dirsum;
    float normsum;
    };


CUDAwqNode *getNode(CUDACTX cuda) {
    CUDAwqNode *work;
    if (cuda->wqFree == NULL) return NULL;
    work = cuda->wqFree;
    cuda->wqFree = work->next;
    work->ctx = cuda;
    work->checkFcn = NULL;
    work->next = cuda->wqCuda;
    work->ppSizeIn = 0;
    work->ppSizeOut = 0;
    work->ppnBuffered = 0;
    work->ppnBlocks = 0;
    return work;
    }

#define NP_ALIGN (128/sizeof(ppResult))
#define NP_ALIGN_MASK (NP_ALIGN-1)


/*
** Occupancy (theoretical):
**   Compute 2.x:  8 (blocks) * PP_THREADS (128) = 1024 (2/3 MAX)
**   Compute 3.x: 16 (blocks) * PP_THREADS (128) = 2048 (MAX!)
**
** To reduce memory, we syncthreads() and flush the particles
** results every SYNC_RATE particles => 8 seems a good choice.
**1
** Shared memory requirements
**  - Particles         32 * nSyncRate (8)                 =  256
**  - warp reduction     4 * nWarps (4) * 32 (threads)     =  512
**  - warp results      24 * nSyncRate (8)  * nWarps (4)   =  768
** TOTAL 1536 * 16 blocks = 24 KB
** Same nSyncRate=16: 2560 * 16 blocks = 40 KB
**
** nvcc -DHAVE_CONFIG_H --ptxas-options=-v -c  -I. -arch=sm_20 cudapp.cu
** ptxas info    : 11 bytes gmem, 8 bytes cmem[14]
** ptxas info    : Compiling entry function '_Z6cudaPPILi4ELi16EEvPK10ppWorkUnitPK7ppInputPK7ILP_BLKP8ppResult' for 'sm_20'
** ptxas info    : Function properties for _Z6cudaPPILi4ELi16EEvPK10ppWorkUnitPK7ppInputPK7ILP_BLKP8ppResult
**     0 bytes stack frame, 0 bytes spill stores, 0 bytes spill loads
** ptxas info    : Used 35 registers, 2560 bytes smem, 64 bytes cmem[0], 12 bytes cmem[16]
**
** Interaction lists limits.
**   Compute 2.x: 65535 * PP_THREADS (128) = 8 million
**   Compute 3.x: effectively unlimited
**
*/

// A good number for nWarps is 4 giving 128 threads per thread block, nSyncRate=8
// Each thread block outputs ay,ay,az,fPot,dirsum,normsum for each particle
template <int nWarps,int nWarpsPerWU,int nSyncRate>
__global__ void cudaInteract(
    const ppWorkUnit * __restrict__ work,
    const ppInput * __restrict__ pPart,
    const ilpBlk<WIDTH> * __restrict__ blk,
    ppResult *out) {
    int i, iSync;
    int iWork, iI, iWarp;

    if (nWarpsPerWU==1) {           // blockDim.z == nWarps, blockDim.y == 1, blockDim.x == 32
        iWork = blockIdx.x * nWarps + threadIdx.z; // Work and corresponding blk
        iI = threadIdx.x; // Index into blk
        iWarp = threadIdx.z;
        }
    else if (nWarps==nWarpsPerWU) { // blockDim.z == 1, blockDim.y == nWarps, blockDim.x == 32
        iWork = blockIdx.x; // Index of work and blk
        iI =   threadIdx.y*blockDim.x + threadIdx.x; // Index of interaction
        iWarp = threadIdx.y;
        }
    else {
        // Calculate our interaction and particle group
        iWork = blockIdx.x*blockDim.z + threadIdx.z; // Work and corresponding blk
        iI =   threadIdx.y*blockDim.x + threadIdx.x; // Thread working on blk
        int iAll = iI + threadIdx.z*blockDim.y*blockDim.x;
        iWarp = iAll / 32;
        }
    int iTinW = iI % 32;

    int nP = work[iWork].nP; // Number of particles
    pPart += work[iWork].iP; // First particle
    int nI = work[iWork].nI; // Number of interactions
    blk += iWork*blockDim.y + threadIdx.y; // blk[threadIdx.x] is our interaction
    out += work[iWork].iO;   // Result for each particle

    __shared__ union {
        ppInput P[nWarps/nWarpsPerWU][nSyncRate];
        float   W[nWarps/nWarpsPerWU][nSyncRate*sizeof(ppInput)/sizeof(float)];
        } Particles;

    __shared__ float reduce[nWarps][32];

    __shared__ float wX[nSyncRate][nWarps];
    __shared__ float wY[nSyncRate][nWarps];
    __shared__ float wZ[nSyncRate][nWarps];

    __shared__ float wPot[nSyncRate][nWarps];
    __shared__ float wDirsum[nSyncRate][nWarps];
    __shared__ float wNormsum[nSyncRate][nWarps];


    // Load the interaction. It is blocked for performance.
    float iX,iY,iZ,iM,ifourh2;
    if (iI < nI) {
        iX = blk->dx[threadIdx.x];
        iY = blk->dy[threadIdx.x];
        iZ = blk->dz[threadIdx.x];
        iM = blk->m[threadIdx.x];
        ifourh2 = blk->fourh2[threadIdx.x];
        }

    // Apply the particles, nSyncRate at a time
    for(iSync=0; iSync<nP; iSync += nSyncRate) {
        int iEnd = nP - iSync;
        if (iEnd > nSyncRate) iEnd=nSyncRate;
        // Preload the bucket of particles - this is a memcpy
        if (iI < iEnd*sizeof(ppInput) / sizeof(float)) {
            Particles.W[threadIdx.z][iI] = (reinterpret_cast<const float *>(pPart+iSync))[iI];
            }
        if (iI < iEnd ) { 
            float ax = Particles.P[threadIdx.z][iI].ax;
            float ay = Particles.P[threadIdx.z][iI].ay;
            float az = Particles.P[threadIdx.z][iI].az;
            Particles.P[threadIdx.z][iI].dImaga = ax*ax + ay*ay + az*az;
            if (Particles.P[threadIdx.z][iI].dImaga > 0.0f)
                Particles.P[threadIdx.z][iI].dImaga = rsqrtf(Particles.P[threadIdx.z][iI].dImaga);
            }
        if (nWarpsPerWU>1) __syncthreads();

        for( i=0; i<iEnd; ++i) {
            float ax=0.0f, ay=0.0f, az=0.0f, fPot=0.0f, dirsum=0.0f, normsum=0.0f;
            if (iI < nI) {
                float fourh2,dir,dir2,dir3;
                float dx = iX + Particles.P[threadIdx.z][i].dx;
                float dy = iY + Particles.P[threadIdx.z][i].dy;
                float dz = iZ + Particles.P[threadIdx.z][i].dz;
                float d2 = dx*dx + dy*dy + dz*dz;
                if (d2 != 0.0f ) { /* Ignore self interactions */
                    fourh2 = ifourh2;
                    if (d2 > fourh2) fourh2 = d2;
                    dir = rsqrtf(fourh2);
                    dir2 = dir*dir;
                    dir3 = dir2*dir;
                    if (d2 < fourh2) {
                        /*
                        ** This uses the Dehnen K1 kernel function now, it's fast!
                        */
                        dir2 *= d2;
                        dir2 = 1.0f - dir2;
                        dir *= 1.0f + dir2*(0.5f + dir2*(3.0f/8.0f + dir2*(45.0f/32.0f)));
                        dir3 *= 1.0f + dir2*(1.5f + dir2*(135.0f/16.0f));
                        }
                    dir3 *= -iM;
                    ax = dx * dir3;
                    ay = dy * dir3;
                    az = dz * dir3;
                    fPot = -iM*dir;
                    /*
                    ** Calculations for determining the timestep.
                    */
                    float adotai;
                    adotai = Particles.P[threadIdx.z][i].ax*ax + Particles.P[threadIdx.z][i].ay*ay + Particles.P[threadIdx.z][i].az*az;
                    if (adotai > 0.0f && d2 >= Particles.P[threadIdx.z][i].fSoft2 ) {
                        adotai *= Particles.P[threadIdx.z][i].dImaga;
                        dirsum = dir*adotai*adotai;
                        normsum = adotai*adotai;
                        }
                    }
                }
            // Horizontal add within each warp -- no sychronization required
            warpReduceAndStore<float,32>(reduce[iWarp],iTinW,ax,           &wX[i][iWarp]);
            warpReduceAndStore<float,32>(reduce[iWarp],iTinW,ay,           &wY[i][iWarp]);
            warpReduceAndStore<float,32>(reduce[iWarp],iTinW,az,           &wZ[i][iWarp]);
            warpReduceAndStore<float,32>(reduce[iWarp],iTinW,fPot,       &wPot[i][iWarp]);
            warpReduceAndStore<float,32>(reduce[iWarp],iTinW,dirsum,  &wDirsum[i][iWarp]);
            warpReduceAndStore<float,32>(reduce[iWarp],iTinW,normsum,&wNormsum[i][iWarp]);
            }

        if (nWarpsPerWU>1) __syncthreads();
        // Assuming four warps & SYNC of 8, the cache looks like this:
        //                0     1     2     3       4
        // Cache: 1x4x8   P0    P0    P0    P0      P1 P1 P1 P1 P2 ... P6 P7 P7 P7 P7
        // Cache: 4x1x8   P0,I0 P0,I1 P0,I2 P0,I3   P1,I0 ... 
        // Cache: 2x2x8   P0,I0 P0,I0 P0,I1 P0,I1   P1,I0
        // every set of 4 threads does another reduce. As long as the
        // number of warps is a power of 2 <= the warp size (32), we
        // can do this step without any further synchronization.
        int nOut = iEnd * nWarpsPerWU; // Normally 64
        if (iI<nOut) {
            int iP    = (iI & ~(nWarpsPerWU-1)) * nWarps/nWarpsPerWU + threadIdx.z; // 0,4,8,...
            int iWarp = iI &  (nWarpsPerWU-1); // 0 .. 3
            int iOut  = iI / nWarpsPerWU + iSync;
         
            warpReduceAndStore<float,nWarpsPerWU>(      &wX[0][0]+iP,iWarp,&out[iOut].ax);
            warpReduceAndStore<float,nWarpsPerWU>(      &wY[0][0]+iP,iWarp,&out[iOut].ay);
            warpReduceAndStore<float,nWarpsPerWU>(      &wZ[0][0]+iP,iWarp,&out[iOut].az);
            warpReduceAndStore<float,nWarpsPerWU>(    &wPot[0][0]+iP,iWarp,&out[iOut].fPot);
            warpReduceAndStore<float,nWarpsPerWU>( &wDirsum[0][0]+iP,iWarp,&out[iOut].dirsum);
            warpReduceAndStore<float,nWarpsPerWU>(&wNormsum[0][0]+iP,iWarp,&out[iOut].normsum);
            }
        }
    }





template <int nWarps,int nWarpsPerWU,int nSyncRate>
__global__ void cudaInteract(
    const ppWorkUnit * __restrict__ work,
    const ppInput * __restrict__ pPart,
    const ilcBlk<WIDTH> * __restrict__ blk,
    ppResult *out) {
    int i, iSync;
    int iWork, iI, iWarp;

    if (nWarpsPerWU==1) {           // blockDim.z == nWarps, blockDim.y == 1, blockDim.x == 32
        iWork = blockIdx.x * nWarps + threadIdx.z; // Work and corresponding blk
        iI = threadIdx.x; // Index into blk
        iWarp = threadIdx.z;
        }
    else if (nWarps==nWarpsPerWU) { // blockDim.z == 1, blockDim.y == nWarps, blockDim.x == 32
        iWork = blockIdx.x; // Index of work and blk
        iI =   threadIdx.y*blockDim.x + threadIdx.x; // Index of interaction
        iWarp = threadIdx.y;
        }
    else {
        // Calculate our interaction and particle group
        iWork = blockIdx.x*blockDim.z + threadIdx.z; // Work and corresponding blk
        iI =   threadIdx.y*blockDim.x + threadIdx.x; // Thread working on blk
        int iAll = iI + threadIdx.z*blockDim.y*blockDim.x;
        iWarp = iAll / 32;
        }
    int iTinW = iI % 32;

    int nP = work[iWork].nP; // Number of particles
    pPart += work[iWork].iP; // First particle
    int nI = work[iWork].nI; // Number of interactions
    blk += iWork*blockDim.y + threadIdx.y; // blk[threadIdx.x] is our interaction
    out += work[iWork].iO;   // Result for each particle

    __shared__ union {
        ppInput P[nWarps/nWarpsPerWU][nSyncRate];
        float   W[nWarps/nWarpsPerWU][nSyncRate*sizeof(ppInput)/sizeof(float)];
        } Particles;

    __shared__ float reduce[nWarps][32];

    __shared__ float wX[nSyncRate][nWarps];
    __shared__ float wY[nSyncRate][nWarps];
    __shared__ float wZ[nSyncRate][nWarps];

    __shared__ float wPot[nSyncRate][nWarps];
    __shared__ float wDirsum[nSyncRate][nWarps];
    __shared__ float wNormsum[nSyncRate][nWarps];


    // Load the interaction. It is blocked for performance.
    float Idx,Idy,Idz;
    float Ixxxx,Ixxxy,Ixxxz,Ixxyz,Ixxyy,Iyyyz,Ixyyz,Ixyyy,Iyyyy;
    float Ixxx,Ixyy,Ixxy,Iyyy,Ixxz,Iyyz,Ixyz;
    float Ixx,Ixy,Ixz,Iyy,Iyz;
#ifdef USE_DIAPOLE
    float Ix,Iy,Iz;
#endif
    float Im,Iu;
    if (iI < nI) {
        Idx = blk->dx[threadIdx.x];
        Idy = blk->dy[threadIdx.x];
        Idz = blk->dz[threadIdx.x];
        Ixxxx = blk->xxxx[threadIdx.x];
        Ixxxy = blk->xxxy[threadIdx.x];
        Ixxxz = blk->xxxz[threadIdx.x];
        Ixxyz = blk->xxyz[threadIdx.x];
        Ixxyy = blk->xxyy[threadIdx.x];
        Iyyyz = blk->yyyz[threadIdx.x];
        Ixyyz = blk->xyyz[threadIdx.x];
        Ixyyy = blk->xyyy[threadIdx.x];
        Iyyyy = blk->yyyy[threadIdx.x];
        Ixxx = blk->xxx[threadIdx.x];
        Ixyy = blk->xyy[threadIdx.x];
        Ixxy = blk->xxy[threadIdx.x];
        Iyyy = blk->yyy[threadIdx.x];
        Ixxz = blk->xxz[threadIdx.x];
        Iyyz = blk->yyz[threadIdx.x];
        Ixyz = blk->xyz[threadIdx.x];
        Ixx = blk->xx[threadIdx.x];
        Ixy = blk->xy[threadIdx.x];
        Ixz = blk->xz[threadIdx.x];
        Iyy = blk->yy[threadIdx.x];
        Iyz = blk->yz[threadIdx.x];
#ifdef USE_DIAPOLE
        Ix = blk->x[threadIdx.x];
        Iy = blk->y[threadIdx.x];
        Iz = blk->z[threadIdx.x];
#endif
        Im = blk->m[threadIdx.x];
        Iu = blk->u[threadIdx.x];
        }
    for(iSync=0; iSync<nP; iSync += nSyncRate) {
        int iEnd = nP - iSync;
        if (iEnd > nSyncRate) iEnd=nSyncRate;
        // Preload the bucket of particles - this is a memcpy
        if (iI < iEnd*sizeof(ppInput) / sizeof(float)) {
            Particles.W[threadIdx.z][iI] = (reinterpret_cast<const float *>(pPart+iSync))[iI];
            }
        if (iI < iEnd ) { 
            float ax = Particles.P[threadIdx.z][iI].ax;
            float ay = Particles.P[threadIdx.z][iI].ay;
            float az = Particles.P[threadIdx.z][iI].az;
            Particles.P[threadIdx.z][iI].dImaga = ax*ax + ay*ay + az*az;
            if (Particles.P[threadIdx.z][iI].dImaga > 0.0f)
                Particles.P[threadIdx.z][iI].dImaga = rsqrtf(Particles.P[threadIdx.z][iI].dImaga);
            }
        if (nWarpsPerWU>1) __syncthreads();

        for( i=0; i<iEnd; ++i) {
            float ax=0.0f, ay=0.0f, az=0.0f, fPot=0.0f, dirsum=0.0f, normsum=0.0f;
            if (iI < nI) {
		const float onethird = 1.0f/3.0f;
                float dx = Idx + Particles.P[threadIdx.z][i].dx;
                float dy = Idy + Particles.P[threadIdx.z][i].dy;
                float dz = Idz + Particles.P[threadIdx.z][i].dz;
		float d2 = dx*dx + dy*dy + dz*dz;
		float dir = rsqrtf(d2);
		float u = Iu*dir;
		float g1 = dir*u;
		float g2 = 3.0f*g1*u;
		float g3 = 5.0f*g2*u;
		float g4 = 7.0f*g3*u;
		/*
		** Calculate the funky distance terms.
		*/
		float x = dx*dir;
		float y = dy*dir;
		float z = dz*dir;
		float xx = 0.5f*x*x;
		float xy = x*y;
		float xz = x*z;
		float yy = 0.5f*y*y;
		float yz = y*z;
		float zz = 0.5f*z*z;
		float xxx = x*(onethird*xx - zz);
		float xxz = z*(xx - onethird*zz);
		float yyy = y*(onethird*yy - zz);
		float yyz = z*(yy - onethird*zz);
		xx -= zz;
		yy -= zz;
		float xxy = y*xx;
		float xyy = x*yy;
		float xyz = xy*z;

	    /*
	    ** Now calculate the interaction up to Hexadecapole order.
	    */
		float tx = g4*(Ixxxx*xxx + Ixyyy*yyy + Ixxxy*xxy + Ixxxz*xxz + Ixxyy*xyy + Ixxyz*xyz + Ixyyz*yyz);
		float ty = g4*(Ixyyy*xyy + Ixxxy*xxx + Iyyyy*yyy + Iyyyz*yyz + Ixxyy*xxy + Ixxyz*xxz + Ixyyz*xyz);
		float tz = g4*(-Ixxxx*xxz - (Ixyyy + Ixxxy)*xyz - Iyyyy*yyz + Ixxxz*xxx + Iyyyz*yyy - Ixxyy*(xxz + yyz) + Ixxyz*xxy + Ixyyz*xyy);
		g4 = 0.25f*(tx*x + ty*y + tz*z);
		xxx = g3*(Ixxx*xx + Ixyy*yy + Ixxy*xy + Ixxz*xz + Ixyz*yz);
		xxy = g3*(Ixyy*xy + Ixxy*xx + Iyyy*yy + Iyyz*yz + Ixyz*xz);
		xxz = g3*(-(Ixxx + Ixyy)*xz - (Ixxy + Iyyy)*yz + Ixxz*xx + Iyyz*yy + Ixyz*xy);
		g3 = onethird*(xxx*x + xxy*y + xxz*z);
		xx = g2*(Ixx*x + Ixy*y + Ixz*z);
		xy = g2*(Iyy*y + Ixy*x + Iyz*z);
		xz = g2*(-(Ixx + Iyy)*z + Ixz*x + Iyz*y);
		g2 = 0.5f*(xx*x + xy*y + xz*z);
		float g0 = dir * Im;
		fPot = -(g0 + g2 + g3 + g4);
		g0 += 5.0f*g2 + 7.0f*g3 + 9.0f*g4;
#ifdef USE_DIAPOLE
		yy = g1*Ix;
		yz = g1*Iy;
		zz = g1*Iz;
		g1 = yy*x + yz*y + zz*z;
		fPot -= g1;
		g0 += 3.0f*g1; 
#else
		yy = 0.0f;
		yz = 0.0f;
		zz = 0.0f;
#endif
		ax = dir*(yy + xx + xxx + tx - x*g0);
		ay = dir*(yz + xy + xxy + ty - y*g0);
		az = dir*(zz + xz + xxz + tz - z*g0);

                /*
                ** Calculations for determining the timestep.
                */
                float adotai;
                adotai = Particles.P[threadIdx.z][i].ax*ax + Particles.P[threadIdx.z][i].ay*ay + Particles.P[threadIdx.z][i].az*az;
                if (adotai > 0.0f && d2 >= Particles.P[threadIdx.z][i].fSoft2 ) {
                    adotai *= Particles.P[threadIdx.z][i].dImaga;
                    dirsum = dir*adotai*adotai;
                    normsum = adotai*adotai;
                    }
                }
            // Horizontal add within each warp -- no sychronization required
            warpReduceAndStore<float,32>(reduce[iWarp],iTinW,ax,     &wX[i][iWarp]);
            warpReduceAndStore<float,32>(reduce[iWarp],iTinW,ay,     &wY[i][iWarp]);
            warpReduceAndStore<float,32>(reduce[iWarp],iTinW,az,     &wZ[i][iWarp]);
            warpReduceAndStore<float,32>(reduce[iWarp],iTinW,fPot,   &wPot[i][iWarp]);
            warpReduceAndStore<float,32>(reduce[iWarp],iTinW,dirsum, &wDirsum[i][iWarp]);
            warpReduceAndStore<float,32>(reduce[iWarp],iTinW,normsum,&wNormsum[i][iWarp]);
	}

        if (nWarpsPerWU>1) __syncthreads();
        // Assuming four warps & SYNC of 8, the cache looks like this:
        //                0     1     2     3       4
        // Cache: 1x4x8   P0    P0    P0    P0      P1 P1 P1 P1 P2 ... P6 P7 P7 P7 P7
        // Cache: 4x1x8   P0,I0 P0,I1 P0,I2 P0,I3   P1,I0 ... 
        // Cache: 2x2x8   P0,I0 P0,I0 P0,I1 P0,I1   P1,I0
        // every set of 4 threads does another reduce. As long as the
        // number of warps is a power of 2 <= the warp size (32), we
        // can do this step without any further synchronization.
        int nOut = iEnd * nWarpsPerWU; // Normally 64
        if (iI<nOut) {
            int iP    = (iI & ~(nWarpsPerWU-1)) * nWarps/nWarpsPerWU + threadIdx.z; // 0,4,8,...
            int iWarp = iI &  (nWarpsPerWU-1); // 0 .. 3
            int iOut  = iI / nWarpsPerWU + iSync;
         
            warpReduceAndStore<float,nWarpsPerWU>(      &wX[0][0]+iP,iWarp,&out[iOut].ax);
            warpReduceAndStore<float,nWarpsPerWU>(      &wY[0][0]+iP,iWarp,&out[iOut].ay);
            warpReduceAndStore<float,nWarpsPerWU>(      &wZ[0][0]+iP,iWarp,&out[iOut].az);
            warpReduceAndStore<float,nWarpsPerWU>(    &wPot[0][0]+iP,iWarp,&out[iOut].fPot);
            warpReduceAndStore<float,nWarpsPerWU>( &wDirsum[0][0]+iP,iWarp,&out[iOut].dirsum);
            warpReduceAndStore<float,nWarpsPerWU>(&wNormsum[0][0]+iP,iWarp,&out[iOut].normsum);
            }
        }
    }

extern "C"
void pkdParticleWorkDone(workParticle *wp);

template<int nIntPerWU>
int CUDAcheckWorkInteraction( void *vpp, void *vwork ) {
    CUDAwqNode *work = reinterpret_cast<CUDAwqNode *>(vwork);
    ppResult *pR       = reinterpret_cast<ppResult *>(work->pHostBuf);
    int ib, iw, ip;

    for( ib=0; ib<work->ppnBuffered; ++ib) {
        workParticle *wp = work->ppWP[ib];
        PINFOOUT *pInfoOut = wp->pInfoOut;
        int nWork = (work->ppNI[ib] + nIntPerWU - 1) / nIntPerWU;
        for(iw=0; iw<nWork; ++iw) {
            for(ip=0; ip<wp->nP; ++ip) {
                pInfoOut[ip].a[0]    += pR[ip].ax;
                pInfoOut[ip].a[1]    += pR[ip].ay;
                pInfoOut[ip].a[2]    += pR[ip].az;
                pInfoOut[ip].fPot    += pR[ip].fPot;
                pInfoOut[ip].dirsum  += pR[ip].dirsum;
                pInfoOut[ip].normsum += pR[ip].normsum;
                }
            pR += wp->nP;
            }
        pkdParticleWorkDone(wp);

        }
    return 0;
    }

extern "C"
void CUDAsetupPP(void) {
    //cudaFuncSetCacheConfig(cudaPP,cudaFuncCachePreferL1);
    }


template<int nIntPerTB, int nIntPerWU, typename BLK>
void CUDA_sendWork(CUDACTX cuda,CUDAwqNode **head) {
    CUDAwqNode *work = *head;
    if (work != NULL) {
        int i, j;
        int iI=0, iP=0, iO=0;
        int nBufferOut = 0;
        int nBlkPer = nIntPerWU / WIDTH;
        int nWork = work->ppnBlocks/nBlkPer;

        // The interation blocks -- already copied to the host memory
        BLK * __restrict__ blkHost = reinterpret_cast<BLK*>(work->pHostBuf);
        BLK * __restrict__ blkCuda = reinterpret_cast<BLK*>(work->pCudaBufIn);
        
        // The interaction block descriptors
        ppWorkUnit * __restrict__ wuHost = reinterpret_cast<ppWorkUnit *>(blkHost + work->ppnBlocks);
        ppWorkUnit * __restrict__ wuCuda = reinterpret_cast<ppWorkUnit *>(blkCuda + work->ppnBlocks);

        // The particle information
        ppInput * __restrict__ partHost = reinterpret_cast<ppInput *>(wuHost + ((nWork+7)&~7));
        ppInput * __restrict__ partCuda = reinterpret_cast<ppInput *>(wuCuda + ((nWork+7)&~7));

        for( i=0; i<work->ppnBuffered; ++i) {
            const int nP = work->ppWP[i]->nP;
            PINFOIN *pInfoIn = work->ppWP[i]->pInfoIn;
            int nPaligned = (nP+NP_ALIGN_MASK) & ~NP_ALIGN_MASK;
            int nInteract = work->ppNI[i];
            int nBlocks = (nInteract+nIntPerWU-1) / nIntPerWU;

            // Generate a interaction block descriptor for each block
            for(j=0; j<nBlocks; ++j) {
                wuHost->nP = nP;
                wuHost->iP = iP;
                wuHost->nI = nInteract > nIntPerWU ? nIntPerWU : nInteract;
                wuHost->iO = iO;
                iO += nP;
                nInteract -= wuHost->nI;
                ++wuHost;
                ++iI;
                nBufferOut += nP * sizeof(ppResult);
                }
            assert(nInteract==0);

            // Copy in nP particles
            for(j=0; j<nP; ++j) {
                partHost[j].dx =  pInfoIn[j].r[0];
                partHost[j].dy =  pInfoIn[j].r[1];
                partHost[j].dz =  pInfoIn[j].r[2];
                partHost[j].ax =  pInfoIn[j].a[0];
                partHost[j].ay =  pInfoIn[j].a[1];
                partHost[j].az =  pInfoIn[j].a[2];
                partHost[j].fSoft2 = pInfoIn[j].fSmooth2;
                /*partHost[j].dImaga = 0;*/
                }
            partHost += nPaligned;
            iP += nPaligned;
            }
        assert(iI == work->ppnBlocks/nBlkPer);
        /* Pad the work out so all work units have valid data */
        int nWUPerTB = nIntPerTB/nIntPerWU;
        while((iI&(nWUPerTB-1)) != 0) {
            wuHost->nP = 0;
            wuHost->iP = 0;
            wuHost->nI = 0;
            wuHost->iO = iO;
            ++wuHost;
            ++iI;
            }
        ppResult *pCudaBufOut = reinterpret_cast<ppResult *>(work->pCudaBufOut);
        int nBufferIn = reinterpret_cast<char *>(partHost) - reinterpret_cast<char *>(work->pHostBuf);
        CUDA_CHECK(cudaMemcpyAsync,(blkCuda, blkHost, nBufferIn, cudaMemcpyHostToDevice, work->stream));

        assert((iI & (nWUPerTB-1)) == 0);
        dim3 dimBlock( WIDTH, nIntPerWU/WIDTH, nIntPerTB/nIntPerWU );
        dim3 dimGrid( iI/nWUPerTB, 1,1);
//        printf("%d x %d x %d X %d x %d x %d\n",
//            dimBlock.x, dimBlock.y, dimBlock.z,
//            dimGrid.x, dimGrid.y, dimGrid.z);
        cudaInteract<WARPS,nIntPerWU/32,SYNC_RATE*nIntPerWU/nIntPerTB>
            <<<dimGrid, dimBlock, 0, work->stream>>>
            (wuCuda,partCuda,blkCuda,pCudaBufOut );

        CUDA_CHECK(cudaMemcpyAsync,(blkHost, work->pCudaBufOut, nBufferOut, cudaMemcpyDeviceToHost, work->stream)
);
        CUDA_CHECK(cudaEventRecord,(work->event,work->stream));

        work->next = cuda->wqCuda;
        cuda->wqCuda = work;

        *head = NULL;
        }
    }

extern "C"
void CUDA_sendWork(void *cudaCtx) {
    CUDACTX cuda = reinterpret_cast<CUDACTX>(cudaCtx);
    CUDA_sendWork< TB_THREADS, PP_WU, ilpBlk<WIDTH> >(cuda,&cuda->nodePP);
    CUDA_sendWork< TB_THREADS, PC_WU, ilcBlk<WIDTH> >(cuda,&cuda->nodePC);
    }


/*
** The following routines copy interactions from an ILP or ILC
** to a CUDA interaction block. The sizes may be different.
*/
template<int n>
int copyBLKs(ilpBlk<n> *out, ILP_BLK *in,int nIlp) {
    assert(n==ILP_PART_PER_BLK);
    int i, nBlk = (nIlp+n-1) / n;
    for(i=0; i<nBlk; ++i) memcpy(&out[i],&in[i],sizeof(out[i]));
    return nBlk;
    }

template<int n>
int copyBLKs(ilcBlk<n> *out, ILC_BLK *in,int nIlp) {
    assert(n==ILC_PART_PER_BLK);
    int i, nBlk = (nIlp+n-1) / n;
    for(i=0; i<nBlk; ++i) memcpy(&out[i],&in[i],sizeof(out[i]));
    return nBlk;
    }

// nIntPer: number of interactions handled per work unit: e.g., 128
template<int nIntPerTB, int nIntPerWU, typename TILE,typename BLK>
int CUDA_queue(CUDACTX cuda,CUDAwqNode **head,workParticle *wp, TILE tile) {
    CUDAwqNode *work = *head;
    /*const int nBlkPerTB = nIntPerTB / WIDTH;*/
    const int nBlkPerWU = nIntPerWU / WIDTH;
    const int nP = wp->nP;
    const int nInteract = tile->lstTile.nBlocks*ILP_PART_PER_BLK + tile->lstTile.nInLast;
    const int nPaligned = (nP+NP_ALIGN_MASK) & ~NP_ALIGN_MASK;
    const int nBlocks = tile->lstTile.nBlocks + (tile->lstTile.nInLast?1:0);
    const int nBlocksAligned = (nBlocks + nBlkPerWU - 1) & ~(nBlkPerWU - 1);
    const int nWork = nBlocksAligned / nBlkPerWU;
    const int nBytesIn = nPaligned * sizeof(ppInput) + nBlocksAligned*sizeof(BLK) + nWork*sizeof(ppWorkUnit);
    const int nBytesOut = nP * sizeof(ppResult) * nWork;

    assert(nWork*nBlkPerWU == nBlocksAligned);

    // Figure out the total amount of space we need, and see if there is enough
    // in the current block. If not, send it to the GPU and get another.
    // Space: nPaligned * 4 (float) * 7 (coordinates)
    //        nBlocks * sizeof(ILP_BLK)
    //
    if (work!=NULL && (work->ppSizeIn + nBytesIn + 8*sizeof(ppWorkUnit) > cuda->inCudaBufSize || work->ppSizeOut + nBytesOut > cuda->outCudaBufSize) ) {
        CUDA_sendWork<nIntPerTB,nIntPerWU,BLK>(cuda,head);
        work = NULL;
        assert(*head==NULL);
        }

    // If we don't have a PP work element, try to grab one
    if (work==NULL) {
        *head = work = getNode(cuda);
        if (work==NULL) return 0;
        work->checkFcn = CUDAcheckWorkInteraction<nIntPerWU>;
        work->ppSizeIn = 0;
        work->ppSizeOut = 0;
        work->ppnBuffered = 0;
        work->ppnBlocks = 0;
        }
    if (work==NULL) return 0;

    // Copy in the interactions. The ILP tiles can then be freed/reused.
    BLK *blk = reinterpret_cast<BLK *>(work->pHostBuf);
    //for(i=0; i<nBlocks; ++i) blk[work->ppnBlocks++] = tile->blk[i];
    copyBLKs(blk+work->ppnBlocks,tile->blk,nInteract);
    work->ppnBlocks += nBlocksAligned;

    work->ppSizeIn += nBytesIn;
    work->ppSizeOut += nBytesOut;
    work->ppNI[work->ppnBuffered] = nInteract;
    work->ppWP[work->ppnBuffered] = wp;
    ++wp->nRefs;

    if ( ++work->ppnBuffered == CUDA_PP_MAX_BUFFERED) CUDA_sendWork<nIntPerTB,nIntPerWU,BLK>(cuda,head);

    return 1;
    }

extern "C"
int CUDA_queuePP(void *cudaCtx,workParticle *wp, ILPTILE tile) {
    CUDACTX cuda = reinterpret_cast<CUDACTX>(cudaCtx);
    return CUDA_queue< TB_THREADS, PP_WU, ILPTILE,ilpBlk<WIDTH> >(cuda,&cuda->nodePP,wp,tile);
    }

extern "C"
int CUDA_queuePC(void *cudaCtx,workParticle *wp, ILCTILE tile) {
    CUDACTX cuda = reinterpret_cast<CUDACTX>(cudaCtx);
    return CUDA_queue< TB_THREADS, PC_WU, ILCTILE,ilcBlk<WIDTH> >(cuda,&cuda->nodePC,wp,tile);
    }

