#ifndef ILC_H
#define ILC_H
#ifdef LOCAL_EXPANSION
#include <stdint.h>

#ifndef ILC_PART_PER_TILE
#define ILC_PART_PER_TILE 1024 /* 1024*100 ~ 100k */
#endif

#if !defined(__CUDACC__)
#include "simd.h"
#endif

/*
** We use a union here so that the compiler can properly align the values.
*/
typedef union {
    float f[ILC_PART_PER_TILE];
#if defined(USE_SIMD_PC) && !defined(__CUDACC__)
    v4sf p[ILC_PART_PER_TILE/SIMD_WIDTH];
#endif
    } ilcFloat;

typedef union {
    uint64_t i[ILC_PART_PER_TILE];
    } ilcInt64;


typedef struct {
    ilcFloat dx,dy,dz;
    ilcFloat xxxx,xxxy,xxxz,xxyz,xxyy,yyyz,xyyz,xyyy,yyyy;
    ilcFloat xxx,xyy,xxy,yyy,xxz,yyz,xyz;
    ilcFloat xx,xy,xz,yy,yz;
    ilcFloat m,u;
    } ILC_BLK;

typedef struct ilcTile {
    struct ilcTile *next;
    struct ilcTile *prev;
    uint32_t nMaxCell;          /* Maximum number of cells in this tile */
    uint32_t nCell;             /* Current number of cells */

  ilcFloat dx,dy,dz;
//  ilcFloat d2;
#ifdef HERMITE
  ilcFloat vx,vy,vz;
#endif
  ilcFloat xxxx,xxxy,xxxz,xxyz,xxyy,yyyz,xyyz,xyyy,yyyy;
  ilcFloat xxx,xyy,xxy,yyy,xxz,yyz,xyz;
  ilcFloat xx,xy,xz,yy,yz;
  ilcFloat m,u;
/*
** ILP has some extra data here that gets partitioned to find the closest
** particles. ILC doesn't need this (at least not yet).
*/
    } *ILCTILE;

typedef struct ilcContext {
    ILCTILE first;              /* first tile in the chain */
    ILCTILE tile;               /* Current tile in the chain */
    double cx, cy, cz;          /* Center coordinates */
    uint32_t nPrevious;         /* Particles in tiles prior to "tile" */

    } *ILC;

typedef struct {
    ILCTILE  tile;
    uint32_t nCell;
    uint32_t nPrevious;
    } ILCCHECKPT;

ILCTILE ilcExtend(ILC ilc);    /* Add tile and return new tile */
ILCTILE ilcClear(ILC ilc);     /* Go back to, and return first tile (empty) */
void ilcInitialize(ILC *ilc);
void ilcFinish(ILC ilc);
size_t ilcMemory(ILC ilc);

static inline void ilcCheckPt(ILC ilc,ILCCHECKPT *cp) {
    cp->tile = ilc->tile;
    cp->nCell = ilc->tile->nCell;
    cp->nPrevious = ilc->nPrevious;
    }

static inline void ilcRestore(ILC ilc,ILCCHECKPT *cp) {
    ilc->tile = cp->tile;
    ilc->nPrevious = cp->nPrevious;
    ilc->tile->nCell = cp->nCell;
    }

static inline uint32_t ilcCount(ILC ilc) {
    return ilc->nPrevious + ilc->tile->nCell;
    }

#if defined(HERMITE)
#define ilcAppend_2(ilc,VX,VY,VZ)					\
    tile->vx.f[ILC_APPEND_i] = (VX);					\
    tile->vy.f[ILC_APPEND_i] = (VY);					\
    tile->vz.f[ILC_APPEND_i] = (VZ);
#else
#define ilcAppend_2(ilc,VX,VY,VZ)
#endif

/*
** The X, Y and Z coordinates must already be relative to cx, cy and cz!
*/
#define ilcAppendFloat(ilc,X,Y,Z,M,U,VX,VY,VZ)				\
    {									\
    ILCTILE tile = (ilc)->tile;						\
    uint_fast32_t ILC_APPEND_i;						\
    if ( tile->nCell == tile->nMaxCell ) tile = ilcExtend((ilc));	\
    ILC_APPEND_i = tile->nCell;						\
    tile->dx.f[ILC_APPEND_i] = (X);			\
    tile->dy.f[ILC_APPEND_i] = (Y);			\
    tile->dz.f[ILC_APPEND_i] = (Z);			\
    assert( (M)->m > 0.0 );						\
    tile->m.f[ILC_APPEND_i] = (M)->m;					\
    tile->u.f[ILC_APPEND_i] = (U);					\
    tile->xx.f[ILC_APPEND_i] = (M)->xx;      				\
    tile->xy.f[ILC_APPEND_i] = (M)->xy;	       			\
    tile->xz.f[ILC_APPEND_i] = (M)->xz;	       			\
    tile->yy.f[ILC_APPEND_i] = (M)->yy;      				\
    tile->yz.f[ILC_APPEND_i] = (M)->yz;				\
    tile->xxx.f[ILC_APPEND_i] = (M)->xxx;            			\
    tile->xxy.f[ILC_APPEND_i] = (M)->xxy;            			\
    tile->xxz.f[ILC_APPEND_i] = (M)->xxz;            			\
    tile->xyy.f[ILC_APPEND_i] = (M)->xyy;            			\
    tile->xyz.f[ILC_APPEND_i] = (M)->xyz;            			\
    tile->yyy.f[ILC_APPEND_i] = (M)->yyy;            			\
    tile->yyz.f[ILC_APPEND_i] = (M)->yyz;            			\
    tile->xxxx.f[ILC_APPEND_i] = (M)->xxxx;                  		\
    tile->xxxy.f[ILC_APPEND_i] = (M)->xxxy;                  		\
    tile->xxxz.f[ILC_APPEND_i] = (M)->xxxz;                  		\
    tile->xxyy.f[ILC_APPEND_i] = (M)->xxyy;                  		\
    tile->xxyz.f[ILC_APPEND_i] = (M)->xxyz;                  		\
    tile->xyyy.f[ILC_APPEND_i] = (M)->xyyy;                  		\
    tile->xyyz.f[ILC_APPEND_i] = (M)->xyyz;                  		\
    tile->yyyy.f[ILC_APPEND_i] = (M)->yyyy;                  		\
    tile->yyyz.f[ILC_APPEND_i] = (M)->yyyz;                  		\
    ilcAppend_2((ilc),VX,VY,VZ);					\
    ++tile->nCell;							\
    }

#define ilcAppend(ilc,X,Y,Z,M,U,VX,VY,VZ)				\
    ilcAppendFloat(ilc,(ilc)->cx-(X),(ilc)->cy-(Y),(ilc)->cz-(Z),M,U,VX,VY,VZ)

#define ILC_LOOP(ilc,ctile) for( ctile=(ilc)->first; ctile!=(ilc)->tile->next; ctile=ctile->next )

#else /* LOCAL_EXPANSION */

#include "moments.h"
/*
** components required for evaluating a multipole interaction.
*/

typedef struct ilCell {
    double x,y,z;
    MOMR mom;
    } ILC;
#endif /* LOCAL_EXPANSION */
#endif
