#ifndef ILP_H
#define ILP_H
#ifdef LOCAL_EXPANSION
#include <stdint.h>
#include <string.h>

#ifndef ILP_PART_PER_TILE
#define ILP_PART_PER_TILE 8192 /* 4096*24 ~ 100k */
#endif
#define ILP_PART_PER_BLK 512
#define ILP_BLK_PER_TILE (ILP_PART_PER_TILE/ILP_PART_PER_BLK)

#if !defined(__CUDACC__)
#include "simd.h"
#endif

/*
** We use a union here so that the compiler can properly align the values.
*/
typedef union {
    float f[ILP_PART_PER_BLK];
#if defined(USE_SIMD_PP) && !defined(__CUDACC__)
    v4sf p[ILP_PART_PER_BLK/SIMD_WIDTH];
#endif
    } ilpVein;

typedef union {
    float f[ILP_PART_PER_TILE];
#if defined(USE_SIMD_PP) && !defined(__CUDACC__)
    v4sf p[ILP_PART_PER_TILE/SIMD_WIDTH];
#endif
    } ilpFloat;

typedef union {
    int64_t i[ILP_PART_PER_TILE]; /* Signed because negative marks softened cells */
    } ilpInt64;


typedef struct {
    ilpVein dx, dy, dz;    /* Offset from ilp->cx, cy, cz */
    ilpVein m;             /* Mass */
    ilpVein fourh2;        /* Softening: calculated */
    } ILP_BLK;

typedef struct ilpTile {
    struct ilpTile *next;
    struct ilpTile *prev;
    uint32_t nMaxPart;          /* Maximum number of particles in this tile */
    uint32_t nPart;             /* Current number of particles */
    ILP_BLK *blk;

/* #ifdef HERMITE */
    ilpFloat vx, vy, vz;
/* #endif */
/* #if defined(SYMBA) || defined(PLANETS) */
    ilpInt64 iOrder;
/* #endif */
    } *ILPTILE;

typedef struct ilpContext {
    ILPTILE first;              /* first tile in the chain */
    ILPTILE tile;               /* Current tile in the chain */
    double cx, cy, cz;          /* Center coordinates */
    uint32_t nPrevious;         /* Particles in tiles prior to "tile" */
#ifdef USE_CACHE
    uint32_t nCache;
    ILP_BLK cache;
#endif
    } *ILP;

typedef struct {
    ILPTILE  tile;
    uint32_t nPart;
    uint32_t nPrevious;
    } ILPCHECKPT;

ILPTILE ilpExtend(ILP ilp);    /* Add tile and return new tile */
ILPTILE ilpClear(ILP ilp);     /* Go back to, and return first tile (empty) */
void ilpInitialize(ILP *ilp);
void ilpFinish(ILP ilp);
size_t ilpMemory(ILP ilp);

static inline void ilpCheckPt(ILP ilp,ILPCHECKPT *cp) {
    cp->tile = ilp->tile;
    cp->nPart = ilp->tile->nPart;
    cp->nPrevious = ilp->nPrevious;
    }

static inline void ilpRestore(ILP ilp,ILPCHECKPT *cp) {
    ilp->tile = cp->tile;
    ilp->nPrevious = cp->nPrevious;
    ilp->tile->nPart = cp->nPart;
    }

static inline uint32_t ilpCount(ILP ilp) {
    return ilp->nPrevious + ilp->tile->nPart;
    }

static inline void ilpFlush(ILP ilp) {
#ifdef USE_CACHE
    if (ilp->nCache) {
	ILPTILE tile = (ilp)->tile;
	uint_fast32_t blk = tile->nPart / ILP_PART_PER_BLK;
	memcpy(tile->blk + blk, &ilp->cache, sizeof(ILP_BLK));
	tile->nPart += ILP_PART_PER_BLK;
	ilp->nCache = 0;
	}
#endif
    }

/*
** The X, Y and Z coordinates must already be relative to cx, cy and cz!
*/

static inline void ilpAppendFloat(ILP ilp, float X, float Y, float Z, float M, float S,
    uint64_t I, float VX, float VY, float VZ ) {
    ILPTILE tile = (ilp)->tile;
    uint_fast32_t i = tile->nPart;
    uint_fast32_t blk = i / ILP_PART_PER_BLK;
    uint_fast32_t prt = i - blk * ILP_PART_PER_BLK;
#ifdef USE_CACHE
    uint_fast32_t ni = ilp->nCache;
    ilp->cache.dx.f[ni] = (X);
    ilp->cache.dy.f[ni] = (Y);
    ilp->cache.dz.f[ni] = (Z);
    ilp->cache.m.f[ni] = (M);
    ilp->cache.fourh2.f[ni] = (S);
#else
    tile->blk[blk].dx.f[prt] = (X);
    tile->blk[blk].dy.f[prt] = (Y);
    tile->blk[blk].dz.f[prt] = (Z);
    tile->blk[blk].m.f[prt] = (M);
    tile->blk[blk].fourh2.f[prt] = (S);
#endif
    assert( (M) > 0.0 );
    tile->iOrder.i[i] = (I);
    tile->vx.f[i] = (VX);
    tile->vy.f[i] = (VY);
    tile->vz.f[i] = (VZ);
#ifdef USE_CACHE
    if ( ++ilp->nCache == ILP_PART_PER_BLK ) {
	uint_fast32_t blk = i / ILP_PART_PER_BLK;
	uint_fast32_t prt = i - blk * ILP_PART_PER_BLK;
	memcpy(tile->blk + blk, &ilp->cache, sizeof(ILP_BLK));
	tile->nPart += ILP_PART_PER_BLK;
	ilp->nCache = 0;
	if ( tile->nPart == tile->nMaxPart ) tile = ilpExtend((ilp));
	}
#else
    if ( ++tile->nPart == tile->nMaxPart ) tile = ilpExtend((ilp));
#endif
    }

static inline void ilpAppend(ILP ilp, double X, double Y, double Z, float M, float S,
    uint64_t I, float VX, float VY, float VZ ) {
    ilpAppendFloat(ilp,(ilp)->cx-(X),(ilp)->cy-(Y),(ilp)->cz-(Z),M,S,I,VX,VY,VZ);
    }
#define ILP_LOOP(ilp,ptile) for( ptile=(ilp)->first; ptile!=(ilp)->tile->next; ptile=ptile->next )

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
