#ifndef CL_H
#define CL_H
#include <stdint.h>

#ifndef CL_PART_PER_TILE
#define CL_PART_PER_TILE 1024 /* 1024*100 ~ 100k */
#endif

#define CL_PART_PER_BLK 1024
#define CL_BLK_PER_TILE (CL_PART_PER_TILE/CL_PART_PER_BLK)

#if !defined(__CUDACC__)
#include "simd.h"
#endif

#include "lst.h"

/*
** We use a union here so that the compiler can properly align the values.
*/
typedef union {
    float f[CL_PART_PER_TILE];
#if defined(USE_SIMD_OPEN) && !defined(__CUDACC__)
    v4sf p[CL_PART_PER_TILE/SIMD_WIDTH];
#endif
    } clFloat;

typedef union {
    int32_t i[CL_PART_PER_TILE];
#if defined(USE_SIMD_OPEN) && !defined(__CUDACC__)
    v4i     p[CL_PART_PER_TILE/SIMD_WIDTH];
    v4sf    pf[CL_PART_PER_TILE/SIMD_WIDTH];
#endif
    } clInt32;

typedef union {
    uint64_t i[CL_PART_PER_TILE];
    } clInt64;

typedef struct {
    clInt32 iOpen;
    clInt32 iCell;
    clInt32 id;
    clInt32 iLower;
    clInt32 nc;
    clFloat cOpen;
    clFloat m;
    clFloat fourh2;
    clFloat x;
    clFloat y;
    clFloat z;
    clFloat xOffset;
    clFloat yOffset;
    clFloat zOffset;
    clFloat xCenter;
    clFloat yCenter;
    clFloat zCenter;
    clFloat xMax;
    clFloat yMax;
    clFloat zMax;
    } CL_BLK;


typedef struct clTile {
    LSTTILE lstTile;
    CL_BLK *blk;
    } *CLTILE;

typedef struct clContext {
    LST lst;
    } *CL;

#define clClear(cl) lstClear(&(cl)->lst)
#define clExtend(cl) lstExtend(&(cl)->lst)
#define clMemory(cl) lstMemory(&(cl)->lst)
#define clCount(cl) lstCount(&(cl)->lst)
#define clClone(dst,src) lstClone(&(dst)->lst,&(src)->lst)

void clInitialize(CL *cl,LSTFREELIST *clFreeList);
void clFinish(CL cl);

static inline void clAppendAll(CL cl,int iCell, int id,int iLower,int nc,
    float cOpen,float m,float fourh2,float x, float y, float z,
    float xOffset,float yOffset,float zOffset,float xCenter,float yCenter,float zCenter,
    float xMax,float yMax,float zMax,int iOpen) {
    CLTILE tile = lstReposition(&cl->lst);
    uint_fast32_t blk = tile->lstTile.nBlocks;
    uint_fast32_t prt = tile->lstTile.nInLast;
    tile->blk[blk].iOpen.i[prt] = iOpen;
    tile->blk[blk].iCell.i[prt] = (iCell);
    tile->blk[blk].id.i[prt] = (id);
    tile->blk[blk].iLower.i[prt] = (iLower);
    tile->blk[blk].nc.i[prt] = (nc);
    tile->blk[blk].cOpen.f[prt] = (cOpen);
    tile->blk[blk].m.f[prt] = (m);
    tile->blk[blk].fourh2.f[prt] = (fourh2);
    tile->blk[blk].x.f[prt] = x;
    tile->blk[blk].y.f[prt] = y;
    tile->blk[blk].z.f[prt] = z;
    tile->blk[blk].xOffset.f[prt] = xOffset;
    tile->blk[blk].yOffset.f[prt] = yOffset;
    tile->blk[blk].zOffset.f[prt] = zOffset;
    tile->blk[blk].xCenter.f[prt] = xCenter;
    tile->blk[blk].yCenter.f[prt] = yCenter;
    tile->blk[blk].zCenter.f[prt] = zCenter;
    tile->blk[blk].xMax.f[prt] = xMax;
    tile->blk[blk].yMax.f[prt] = yMax;
    tile->blk[blk].zMax.f[prt] = zMax;
    ++tile->lstTile.nInLast;
    }

#define clAppend(cl,iCell,id,iLower,nc,cOpen,m,fourh2,r,fOffset,fCenter,fMax)\
    clAppendAll(cl,iCell,id,iLower,nc,cOpen,m,fourh2,r[0],r[1],r[2],\
    fOffset[0],fOffset[1],fOffset[2],fCenter[0],fCenter[1],fCenter[2],\
    fMax[0],fMax[1],fMax[2],0)

static inline void clAppendItem(CL cl, CL_BLK *B, int Bi) {
    clAppendAll(cl,B->iCell.i[Bi],B->id.i[Bi],B->iLower.i[Bi], B->nc.i[Bi], B->cOpen.f[Bi], B->m.f[Bi], B->fourh2.f[Bi],
	B->x.f[Bi], B->y.f[Bi], B->z.f[Bi],B->xOffset.f[Bi],B->yOffset.f[Bi],B->zOffset.f[Bi],
	B->xCenter.f[Bi],B->yCenter.f[Bi],B->zCenter.f[Bi],B->xMax.f[Bi],B->yMax.f[Bi],B->zMax.f[Bi],B->iOpen.i[Bi]);
    }

#define CL_LOOP(CL,CL_TILE) for( CL_TILE=(CLTILE)((CL)->lst.list); CL_TILE!=NULL; CL_TILE=(CLTILE)(CL_TILE->lstTile.next))

#endif
