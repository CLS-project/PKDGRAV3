#ifndef ILP_H
#define ILP_H
#ifdef LOCAL_EXPANSION
#include <stdint.h>
#include <string.h>

#include "lst.h"

#define ILP_TILE_SIZE (100*1024) /* 100k */
#ifndef ILP_PART_PER_BLK
#define ILP_PART_PER_BLK (64)
#endif

#if !defined(__CUDACC__)
#include "simd.h"
#endif

/*
** We use a union here so that the compiler can properly align the values.
*/
typedef union {
    float f[ILP_PART_PER_BLK];
#if defined(USE_SIMD_PP) && !defined(__CUDACC__)
    v_sf p[ILP_PART_PER_BLK/SIMD_WIDTH];
#endif
    } ilpFloat;

typedef union {
    int64_t i[ILP_PART_PER_BLK]; /* Signed because negative marks softened cells */
    } ilpInt64;

typedef struct {
    ilpFloat dx, dy, dz;    /* Offset from ilp->cx, cy, cz */
    ilpFloat m;             /* Mass */
    ilpFloat fourh2;        /* Softening: calculated */
    } ILP_BLK;

typedef struct {
    ilpFloat vx, vy, vz;
    ilpInt64 iOrder;
    } ILP_EXTRA;

typedef struct ilpTile {
    LSTTILE lstTile;
    ILP_BLK *blk;
    ILP_EXTRA *xtr;
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
    tile->xtr[blk].iOrder.i[prt] = (I);
    tile->xtr[blk].vx.f[prt] = (VX);
    tile->xtr[blk].vy.f[prt] = (VY);
    tile->xtr[blk].vz.f[prt] = (VZ);
    ++tile->lstTile.nInLast;
    }

static inline void ilpAppend(ILP ilp, double X, double Y, double Z, float M, float S,
    uint64_t I, float VX, float VY, float VZ ) {
    ilpAppendFloat(ilp,(ilp)->cx-(X),(ilp)->cy-(Y),(ilp)->cz-(Z),M,S,I,VX,VY,VZ);
    }
#define ILP_LOOP(ilp,ptile) for( ptile=(ILPTILE)((ilp)->lst.list); ptile!=NULL; ptile=(ILPTILE)(ptile->lstTile.next) )

#else /* LOCAL_EXPANSION */

typedef struct ilPart {
    double m,x,y,z;
#if defined(SOFTLINEAR)
    double h;
#elif defined(SOFTSQUARE)
    double twoh2;
#else
    double fourh2;
#endif  
#if defined(SYMBA) || defined(PLANETS) || !defined(LOCAL_EXPANSION)
    uint64_t iOrder;
#endif
#if defined(HERMITE) || !defined(LOCAL_EXPANSION)
    double vx,vy,vz;
#endif
    } ILP;
#endif /* LOCAL_EXPANSION */

#endif
