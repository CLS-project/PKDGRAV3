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

#ifndef ILP_H
#define ILP_H
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>

#include "lst.h"

#define ILP_TILE_SIZE (100*1024) /* 100k */
#ifndef ILP_PART_PER_BLK
#define ILP_PART_PER_BLK (32) /* Don't mess with this: see CUDA */
#endif

#if !defined(__CUDACC__)
#include "core/simd.h"
#endif

/*
** We use a union here so that the compiler can properly align the values.
*/
typedef union {
    float f[ILP_PART_PER_BLK];
#if !defined(__CUDACC__)
    v_sf p[ILP_PART_PER_BLK/SIMD_WIDTH];
#endif
    } ilpFloat;

typedef union {
    double d[ILP_PART_PER_BLK];
    } ilpDouble;

typedef union {
    int64_t i[ILP_PART_PER_BLK]; /* Signed because negative marks softened cells */
    } ilpInt64;

typedef struct {
    ilpFloat dx, dy, dz;    /* Offset from ilp->cx, cy, cz */
    ilpFloat m;             /* Mass */
    ilpFloat fourh2;        /* Softening: calculated */
    } ILP_BLK;

#ifdef TIMESTEP_CRITICAL
typedef struct {
    ilpDouble vx, vy, vz;
    ilpInt64 iOrder;
    } ILP_EXTRA;
#endif

typedef struct ilpTile {
    LSTTILE lstTile;
    ILP_BLK *blk;
#ifdef TIMESTEP_CRITICAL
    ILP_EXTRA *xtr;
#endif
    } *ILPTILE;

typedef struct ilpContext {
    LST lst;
    double cx, cy, cz;          /* Center coordinates */
    } *ILP;

#define ILPCHECKPT LSTCHECKPT

void ilpInitialize(ILP *ilp);
void ilpFinish(ILP ilp);

#define ilpClear(ilp) lstClear(&(ilp)->lst)
#define ilpExtend(ilp) lstExtend(&(ilp)->lst)
#define ilpMemory(ilp) lstMemory(&(ilp)->lst)
#define ilpCheckPt(ilp,cp) lstCheckPt(&(ilp)->lst,(cp))
#define ilpRestore(ilp,cp) lstRestore(&(ilp)->lst,(cp))
#define ilpCount(ilp) lstCount(&(ilp)->lst)

/*
** The X, Y and Z coordinates must already be relative to cx, cy and cz!
*/

static inline void ilpAppendFloat(ILP ilp, float X, float Y, float Z, float M, float S,
    uint64_t I, float VX, float VY, float VZ ) {
    ILPTILE tile = (ILPTILE)lstReposition(&ilp->lst);
    uint_fast32_t blk = tile->lstTile.nBlocks;
    uint_fast32_t prt = tile->lstTile.nInLast;
    tile->blk[blk].dx.f[prt] = (X);
    tile->blk[blk].dy.f[prt] = (Y);
    tile->blk[blk].dz.f[prt] = (Z);
    tile->blk[blk].m.f[prt] = (M);
    tile->blk[blk].fourh2.f[prt] = (S);
    assert( (M) > 0.0 );
#ifdef TIMESTEP_CRITICAL
    tile->xtr[blk].iOrder.i[prt] = (I);
    tile->xtr[blk].vx.d[prt] = (VX);
    tile->xtr[blk].vy.d[prt] = (VY);
    tile->xtr[blk].vz.d[prt] = (VZ);
#endif
    ++tile->lstTile.nInLast;
    }

static inline void ilpAppend(ILP ilp, double X, double Y, double Z, float M, float S,
    uint64_t I, float VX, float VY, float VZ ) {
    ilpAppendFloat(ilp,(float)((ilp)->cx-(X)),(float)((ilp)->cy-(Y)),(float)((ilp)->cz-(Z)),M,S,I,VX,VY,VZ);
    }
#define ILP_LOOP(ilp,ptile) for( ptile=(ILPTILE)((ilp)->lst.list); ptile!=NULL; ptile=(ILPTILE)(ptile->lstTile.next) )

#endif
