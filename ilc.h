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

#ifndef ILC_H
#define ILC_H
#include <stdint.h>
#include <assert.h>
#include <stdint.h>

#include "lst.h"

#define ILC_TILE_SIZE (100*1024) /* 100k */
#ifndef ILC_PART_PER_BLK
#define ILC_PART_PER_BLK (32) /* Don't mess with this: see CUDA */
#endif

#if !defined(__CUDACC__)
#include "core/simd.h"
#endif

#include "gravity/moments.h"

/*
** We use a union here so that the compiler can properly align the values.
*/
typedef union {
    float f[ILC_PART_PER_BLK];
#if !defined(__CUDACC__)
    v_sf p[ILC_PART_PER_BLK / SIMD_WIDTH];
#endif
    } ilcFloat;

typedef struct {
    ilcFloat dx,dy,dz;
    ilcFloat xxxx,xxxy,xxxz,xxyz,xxyy,yyyz,xyyz,xyyy,yyyy;
    ilcFloat xxx,xyy,xxy,yyy,xxz,yyz,xyz;
    ilcFloat xx,xy,xz,yy,yz;
    ilcFloat x,y,z;
    ilcFloat m,u;
    } ILC_BLK;

typedef struct ilcTile {
    LSTTILE lstTile;
    ILC_BLK *blk;
    } *ILCTILE;

typedef struct ilcContext {
    LST lst;
    double cx, cy, cz;          /* Center coordinates */
    } *ILC;

#define ILCCHECKPT LSTCHECKPT

#define ilcClear(ilc) lstClear(&(ilc)->lst)
#define ilcExtend(ilc) lstExtend(&(ilc)->lst)
#define ilcMemory(ilc) lstMemory(&(ilc)->lst)
#define ilcCheckPt(ilc,cp) lstCheckPt(&(ilc)->lst,(cp))
#define ilcRestore(ilc,cp) lstRestore(&(ilc)->lst,(cp))
#define ilcCount(ilc) lstCount(&(ilc)->lst)

void ilcInitialize(ILC *ilc);
void ilcFinish(ILC ilc);

/*
** The X, Y and Z coordinates must already be relative to cx, cy and cz!
*/
static inline void ilcAppendFloat(ILC ilc,float X,float Y,float Z,FMOMR *M,float U) {
    ILCTILE tile = (ILCTILE)lstReposition(&ilc->lst);
    uint_fast32_t blk = tile->lstTile.nBlocks;
    uint_fast32_t prt = tile->lstTile.nInLast;

    tile->blk[blk].dx.f[prt] = (X);
    tile->blk[blk].dy.f[prt] = (Y);
    tile->blk[blk].dz.f[prt] = (Z);
    assert( (M)->m > 0.0 );
    tile->blk[blk].m.f[prt] = (M)->m;
    tile->blk[blk].u.f[prt] = (U);
    tile->blk[blk].x.f[prt] = (M)->x;
    tile->blk[blk].y.f[prt] = (M)->y;
    tile->blk[blk].z.f[prt] = (M)->z;
    tile->blk[blk].xx.f[prt] = (M)->xx;
    tile->blk[blk].xy.f[prt] = (M)->xy;
    tile->blk[blk].xz.f[prt] = (M)->xz;
    tile->blk[blk].yy.f[prt] = (M)->yy;
    tile->blk[blk].yz.f[prt] = (M)->yz;
    tile->blk[blk].xxx.f[prt] = (M)->xxx;
    tile->blk[blk].xxy.f[prt] = (M)->xxy;
    tile->blk[blk].xxz.f[prt] = (M)->xxz;
    tile->blk[blk].xyy.f[prt] = (M)->xyy;
    tile->blk[blk].xyz.f[prt] = (M)->xyz;
    tile->blk[blk].yyy.f[prt] = (M)->yyy;
    tile->blk[blk].yyz.f[prt] = (M)->yyz;
    tile->blk[blk].xxxx.f[prt] = (M)->xxxx;
    tile->blk[blk].xxxy.f[prt] = (M)->xxxy;
    tile->blk[blk].xxxz.f[prt] = (M)->xxxz;
    tile->blk[blk].xxyy.f[prt] = (M)->xxyy;
    tile->blk[blk].xxyz.f[prt] = (M)->xxyz;
    tile->blk[blk].xyyy.f[prt] = (M)->xyyy;
    tile->blk[blk].xyyz.f[prt] = (M)->xyyz;
    tile->blk[blk].yyyy.f[prt] = (M)->yyyy;
    tile->blk[blk].yyyz.f[prt] = (M)->yyyz;

    ++tile->lstTile.nInLast;
    }

static inline void ilcAppend(ILC ilc,double X,double Y,double Z,FMOMR *M,float U) {
    ilcAppendFloat(ilc,(float)((ilc)->cx-(X)),(float)((ilc)->cy-(Y)),(float)((ilc)->cz-(Z)),M,U);
    }

#define ILC_LOOP(ilc,ctile) for( ctile=(ILCTILE)((ilc)->lst.list); ctile!=NULL; ctile=(ILCTILE)(ctile->lstTile.next) )

#endif
