/*  This file is part of PKDGRAV3 (http://www.pkdgrav.org/).
 *  Copyright (c) 2001-2018 Joachim Stadel & Douglas Potter
 *
 *  PKDGRAV3 is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  PKDGRAV3 is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with PKDGRAV3.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef CL_H
#define CL_H
#include <stdint.h>

#ifndef CL_PART_PER_TILE
    #define CL_PART_PER_TILE 1024 /* 1024*100 ~ 100k */
#endif

#define CL_PART_PER_BLK 16
#define CL_BLK_PER_TILE (CL_PART_PER_TILE/CL_PART_PER_BLK)

#if !defined(__CUDACC__)
    #include "core/simd.h"
#endif

#include "lst.h"
#include "../SPHOptions.h"

/*
** We use a union here so that the compiler can properly align the values.
*/
typedef union {
    float f[CL_PART_PER_BLK];
#if !defined(__CUDACC__)
    v_sf p[CL_PART_PER_BLK/SIMD_WIDTH];
#endif
} clFloat;

typedef union {
    int32_t i[CL_PART_PER_BLK];
#if !defined(__CUDACC__)
    v_i     p[CL_PART_PER_BLK/SIMD_WIDTH];
#endif
} clInt32;

typedef union {
    uint64_t i[CL_PART_PER_BLK];
} clInt64;

typedef struct {
    clInt32 iOpen;
    clInt32 iCache;
    clInt32 idCell;
    clInt32 iCell;
    clInt32 idLower;
    clInt32 iLower;
    clInt32 idUpper;
    clInt32 iUpper;
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
#if SPHBALLOFBALLS
    clFloat fBoBr2;
    clFloat fBoBxCenter;
    clFloat fBoByCenter;
    clFloat fBoBzCenter;
#endif
#if SPHBOXOFBALLS
    clFloat fBoBxMin;
    clFloat fBoBxMax;
    clFloat fBoByMin;
    clFloat fBoByMax;
    clFloat fBoBzMin;
    clFloat fBoBzMax;
#endif
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

#ifdef __cplusplus
extern "C" {
#endif
void clInitialize(CL *cl,LSTFREELIST *clFreeList);
void clDestroy(CL cl);
#ifdef __cplusplus
}
#endif

static inline void clAppendAll(
    CL cl,int iCache,int idCell,int iCell, int idLower,int iLower,int idUpper,int iUpper,int nc,
    float cOpen,float m,float fourh2,float x, float y, float z,
    float xOffset,float yOffset,float zOffset,float xCenter,float yCenter,float zCenter,
#if SPHBALLOFBALLS
    float xMax,float yMax,float zMax,int iOpen,float fBoBr2,float fBoBxCenter,float fBoByCenter,float fBoBzCenter) {
#endif
#if SPHBOXOFBALLS
    float xMax,float yMax,float zMax,int iOpen,float fBoBxMin,float fBoBxMax,float fBoByMin,float fBoByMax,float fBoBzMin,float fBoBzMax) {
#endif
        CLTILE tile = (CLTILE)lstReposition(&cl->lst);
        uint_fast32_t blk = tile->lstTile.nBlocks;
        uint_fast32_t prt = tile->lstTile.nInLast;
        tile->blk[blk].iCache.i[prt] = iCache;
        tile->blk[blk].iOpen.i[prt] = iOpen;
        tile->blk[blk].idCell.i[prt] = (idCell);
        tile->blk[blk].iCell.i[prt] = (iCell);
        tile->blk[blk].idLower.i[prt] = (idLower);
        tile->blk[blk].iLower.i[prt] = (iLower);
        tile->blk[blk].idUpper.i[prt] = (idUpper);
        tile->blk[blk].iUpper.i[prt] = (iUpper);
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
#if SPHBALLOFBALLS
        tile->blk[blk].fBoBr2.f[prt] = fBoBr2;
        tile->blk[blk].fBoBxCenter.f[prt] = fBoBxCenter;
        tile->blk[blk].fBoByCenter.f[prt] = fBoByCenter;
        tile->blk[blk].fBoBzCenter.f[prt] = fBoBzCenter;
#endif
#if SPHBOXOFBALLS
        tile->blk[blk].fBoBxMin.f[prt] = fBoBxMin;
        tile->blk[blk].fBoBxMax.f[prt] = fBoBxMax;
        tile->blk[blk].fBoByMin.f[prt] = fBoByMin;
        tile->blk[blk].fBoByMax.f[prt] = fBoByMax;
        tile->blk[blk].fBoBzMin.f[prt] = fBoBzMin;
        tile->blk[blk].fBoBzMax.f[prt] = fBoBzMax;
#endif
        ++tile->lstTile.nInLast;
    }

#if SPHBALLOFBALLS
#define clAppend(cl,iCache,idCell,iCell,idLower,iLower,idUpper,iUpper,nc,cOpen,m,fourh2,r,fOffset,fCenter,fMax,fBoBr2,fBoBxCenter,fBoByCenter,fBoBzCenter) \
    clAppendAll(cl,iCache,idCell,iCell,idLower,iLower,idUpper,iUpper,nc,cOpen,m,fourh2,r[0],r[1],r[2], \
    fOffset[0],fOffset[1],fOffset[2],fCenter[0],fCenter[1],fCenter[2],\
    fMax[0],fMax[1],fMax[2],0,fBoBr2,fBoBxCenter,fBoByCenter,fBoBzCenter)
#endif
#if SPHBOXOFBALLS
#define clAppend(cl,iCache,idCell,iCell,idLower,iLower,idUpper,iUpper,nc,cOpen,m,fourh2,r,fOffset,fCenter,fMax,fBoBxMin,fBoBxMax,fBoByMin,fBoByMax,fBoBzMin,fBoBzMax) \
    clAppendAll(cl,iCache,idCell,iCell,idLower,iLower,idUpper,iUpper,nc,cOpen,m,fourh2,r[0],r[1],r[2], \
    fOffset[0],fOffset[1],fOffset[2],fCenter[0],fCenter[1],fCenter[2],\
    fMax[0],fMax[1],fMax[2],0,fBoBxMin,fBoBxMax,fBoByMin,fBoByMax,fBoBzMin,fBoBzMax)
#endif


    static inline void clAppendItem(CL cl, CL_BLK *B, int Bi) {
        clAppendAll(cl,B->iCache.i[Bi],B->idCell.i[Bi],B->iCell.i[Bi],
                    B->idLower.i[Bi],B->iLower.i[Bi],B->idUpper.i[Bi],B->iUpper.i[Bi],
                    B->nc.i[Bi], B->cOpen.f[Bi], B->m.f[Bi], B->fourh2.f[Bi],
                    B->x.f[Bi], B->y.f[Bi], B->z.f[Bi],B->xOffset.f[Bi],B->yOffset.f[Bi],B->zOffset.f[Bi],
#if SPHBALLOFBALLS
                    B->xCenter.f[Bi],B->yCenter.f[Bi],B->zCenter.f[Bi],B->xMax.f[Bi],B->yMax.f[Bi],B->zMax.f[Bi],B->iOpen.i[Bi],B->fBoBr2.f[Bi],B->fBoBxCenter.f[Bi],B->fBoByCenter.f[Bi],B->fBoBzCenter.f[Bi]);
#endif
#if SPHBOXOFBALLS
        B->xCenter.f[Bi],B->yCenter.f[Bi],B->zCenter.f[Bi],B->xMax.f[Bi],B->yMax.f[Bi],B->zMax.f[Bi],B->iOpen.i[Bi],B->fBoBxMin.f[Bi],B->fBoBxMax.f[Bi],B->fBoByMin.f[Bi],B->fBoByMax.f[Bi],B->fBoBzMin.f[Bi],B->fBoBzMax.f[Bi]);
#endif
    }

#define CL_LOOP(CL,CL_TILE) for( CL_TILE=(CLTILE)((CL)->lst.list); CL_TILE!=NULL; CL_TILE=(CLTILE)(CL_TILE->lstTile.next))

#endif
